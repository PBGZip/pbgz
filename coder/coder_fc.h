#include "coder_io.h"
#include "coder.h"
#include "../manager.h"
#include "fc/pp.h"
#include "fc/transform.h"
#include "fc/fc.h"

#define FC_MIN_LEN 32
// #define FC_MIN_LEN 2048 // 先设定保守点的值，因为如果太小时编码后长度可能大于原始长度，这样导致fc压缩报错上面流程不好处理
#define FC_MAX_LEN 2146435072

class coder_fc : public coder
{
public:
    coder_fc(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_fc");
        flushed = false;
        lzp_valid = true;

        // if (this->io->get_level() == 1)
        if (true)
        {
            this->ppHashSize = 17;
            this->ppMinLen = 28;
            // 29 -> 123616168, 比28好
            // 27 -> 75368634
            // 26 -> 75365398
            // 20 -> 75350908
            // 18 -> 75347944
            // 17 -> 75361248
            // 16 -> 75540925
        }
        else
            check_exit(false, ERR_INTERNEL, "Currently level %d is not enabled", this->io->get_level());

        check_exit(!(ppMinLen < 4 || ppMinLen > 255), ERR_INTERNEL, "check failed (%d) : %d", __LINE__, ppMinLen);
        check_exit(!(ppHashSize < 10 || ppHashSize > 28), ERR_INTERNEL, "check failed (%d) : %d", __LINE__, ppMinLen);
    }

    virtual ~coder_fc()
    {
        if (io->m == coder_io::MENC && !flushed)
            encode_flush();
    }

    void encode_line(const uint8_t *in, const int32_t in_len, bool need2hold = false)
    {
        check_exit(io->m != coder_io::MENC, ERR_INTERNEL, "only support block compress, not support line method"); // 暂未扩展成line方式
        check_exit(in_len > FC_MIN_LEN && in_len < FC_MAX_LEN, ERR_INTERNEL, "check failed (%d) : %d", __LINE__, in_len);

        io->m = coder_io::MENC;

        // fprintf(stderr, "preprocess...");

        // 先做pp prehandler
        uint8_t *lout;
        safe_alloc(in_len, uint8_t, lout);
        int lzSize;
        int result = fc_preprocess(in, in + in_len, lout + 1, lout + in_len - 1, ppHashSize, ppMinLen);
        if (result >= 0)
        {
            result = (lout[0] = 1, result + 1);
            lzSize = result;
        }
        else {
            memcpy(lout, in, in_len);
            lzSize = in_len;
            lzp_valid = false;
        }
        check_exit(lzSize > 0, ERR_INTERNEL, "check failed (%d) : %d, in_len %d", __LINE__, lzSize, in_len);

        // 再做bwt
        index = fc_transform(lout, lzSize, &num_indexes, indexes);
        if (in_len < 64 * 1024)
            num_indexes = 0;
        check_exit(index >= FC_OK, ERR_INTERNEL, "check failed (%d) : %d, in_len %d", __LINE__, index, in_len);

        // 再压缩
        uint8_t *fc_buff;
        bool need_alloc = (lzSize + 4096) > in_len;
        if (need_alloc)
        {
            safe_alloc(lzSize + 4096, uint8_t, fc_buff);
        }
        else
            fc_buff = this->io->data;

        // fprintf(stderr, "done\n");
        // fprintf(stderr, "compress...");
        result = fcinit();
        check_exit(FC_OK == result, ERR_INTERNEL, "check failed (%d) : %d, in_len %d", __LINE__, result, in_len);
        result = fc_encode(lout, fc_buff + 1, lzSize, lzSize - 1);
        check_exit(result >= 0, ERR_INTERNEL, "check failed (%d) : %d, in_len %d", __LINE__, result, in_len);
        result = (fc_buff[0] = 1, result + 1);

        if (need_alloc)
        {
            memcpy(this->io->data, fc_buff, result);
            free(fc_buff);
        }
        this->io->data_len += result;
        free(lout);
        // fprintf(stderr, "done\n");
    }

    void encode_flush()
    {
        if (io->m != coder_io::MENC || flushed)
            return;

        io->meta["magic"] = io->get_magic();
        io->meta["bi"] = (Json::Value::Int)(index);
        io->meta["bn"] = (Json::Value::Int)(num_indexes);
        io->meta["lv"] = (Json::Value::Int)((lzp_valid ? 1 : 0));

        memcpy(io->data + this->io->data_len, indexes, 4 * num_indexes);
        this->io->data_len += 4 * num_indexes;

        flushed = true;
    }

    /* 对外解压接口，返回实际解压出来的长度，当解压遇到split_ch时退出 */
    int32_t decode_line(uint8_t *out, int32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold = false)
    {
        check_exit(io->m != coder_io::MDEC, ERR_INTERNEL, "only support block decompress, not support line method"); // 暂未扩展成line方式
        check_exit(split_ch == UINT8_MAX, ERR_INTERNEL, "check failed (%d) : %d", __LINE__, split_ch);

        // fprintf(stderr, "decompress...");
        int result = fcinit();
        check_exit(FC_OK == result, ERR_INTERNEL, "check failed (%d) : %d", __LINE__, result);

        uint8_t *lout;
        safe_alloc(out_len, uint8_t, lout);

        result = fc_decode(io->data + 1, lout);
        int lzSize = result;
        check_exit(lzSize > 0 && lzSize <= out_len, ERR_INTERNEL, "check failed (%d) : %d", __LINE__, lzSize);

        index = io->meta["coder"]["bi"].asInt();
        num_indexes = io->meta["coder"]["bn"].asInt();
        lzp_valid = io->meta["coder"]["lv"].asInt();
        memcpy(indexes, io->data + io->meta["tot_dstlen"].asInt() - (num_indexes * 4), (num_indexes * 4));

        // fprintf(stderr, "done\n");
        // fprintf(stderr, "untransform...");
        result = fc_untransform(lout, lzSize, index, num_indexes, indexes);
        check_exit(result >= FC_OK, ERR_INTERNEL, "check failed (%d) : %d", __LINE__, result);

        if (lzp_valid)
        {
            result = fc_unpreprocess(lout + 1, lout + lzSize, out, ppHashSize, ppMinLen);
            check_exit(result >= 0, ERR_INTERNEL, "check failed (%d) : %d", __LINE__, result);
        } else {
            memcpy(out, lout, lzSize);
            result = lzSize;
        }

        if (lout)
            free(lout);
        io->m = coder_io::MDEC;
        // fprintf(stderr, "done\n");

        return result;
    }

private:
    int ppHashSize;
    int ppMinLen;
    bool flushed;

    // used lzp
    bool lzp_valid;

    // bwt
    int index = FC_BAD_ARGS;
    uint8_t num_indexes;
    int indexes[256];
};