#include "coder_json.h"
#include "../manager.h"

coder_json::coder_json(int32_t level) : cLevel(level)
{
}

coder_json::~coder_json()
{
}

/* 压缩json数据，返回压缩后的长度，如果为负数说明out空间不够，如返回-28,说明out_len长度还需再加28
 * 调该函数则假设out_len肯定是足够的，否则一旦out_len不足则会报错退出
 */
int64_t coder_json::encoder(const Json::Value &in, uint8_t *out, const int64_t out_len)
{
    Json::StreamWriterBuilder builder;
    if (in.empty())
        return 0;
    const std::string encjson = Json::writeString(builder, in);
    ZSTD_CCtx *const cctx = ZSTD_createCCtx();
    check_exit(cctx != NULL, ERR_INTERNEL, "ZSTD_createCCtx() failed!");

    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, cLevel), ERR_INTERNEL, "ZSTD_CCtx_setParameter cLevel failed!");
    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1), ERR_INTERNEL, "ZSTD_CCtx_setParameter checksum flag failed!");

    ZSTD_EndDirective const mode = ZSTD_e_end;
    ZSTD_inBuffer input = {encjson.c_str(), encjson.length(), 0};
    ZSTD_outBuffer output = {out, (size_t)(out_len), 0};
    size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    if (remaining != 0)
    {
        ZSTD_freeCCtx(cctx);
        fprintf(stderr, "ZSTD_compressStream2 error : out buffer lack of %ld bytes\n", remaining);
        manage::instance().exit(ERR_INTERNEL);
    }
    if (input.pos != input.size)
    {
        fprintf(stderr, "ZSTD_compressStream2 error : input.pos %ld != input.size %ld\n", input.pos, input.size);
        manage::instance().exit(ERR_INTERNEL);
    }

    ZSTD_freeCCtx(cctx);
    return output.pos;
}

/* 压缩json数据，在内部申请空间，压缩后数据存在out中 */
void coder_json::encoder(const Json::Value &in, std::string &out)
{
    out.clear();
    Json::StreamWriterBuilder builder;
    if (in.empty())
        return;
    const std::string encjson = Json::writeString(builder, in);
    ZSTD_CCtx *const cctx = ZSTD_createCCtx();
    check_exit(cctx != NULL, ERR_INTERNEL, "ZSTD_createCCtx() failed!");

    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, cLevel), ERR_INTERNEL, "ZSTD_CCtx_setParameter cLevel failed!");
    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1), ERR_INTERNEL, "ZSTD_CCtx_setParameter checksum flag failed!");

    ZSTD_EndDirective const mode = ZSTD_e_end;
    ZSTD_inBuffer input = {encjson.c_str(), encjson.length(), 0};
    uint8_t *buffOut;
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    safe_alloc(buffOutSize, uint8_t, buffOut);

    int32_t finished;
    do
    {
        /* Compress into the output buffer and write all of the output to
             * the file so we can reuse the buffer next iteration.
             */
        ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
        size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
        if (ZSTD_isError(remaining))
        {
            free(buffOut);
            fprintf(stderr, "ZSTD_compressStream2 error : %s \n",
                    ZSTD_getErrorName(remaining));
            manage::instance().exit(ERR_INTERNEL);
        }
        out += std::string((char *)buffOut, output.pos);
        /* If we're on the last chunk we're finished when zstd returns 0,
             * which means its consumed all the input AND finished the frame.
             * Otherwise, we're finished when we've consumed all the input.
             */
        finished = (remaining == 0);
    } while (!finished);

    if (input.pos != input.size)
    {
        free(buffOut);
        fprintf(stderr, "ZSTD_compressStream2 error : input.pos %ld != input.size %ld\n", input.pos, input.size);
        manage::instance().exit(ERR_INTERNEL);
    }

    ZSTD_freeCCtx(cctx);
    free(buffOut);
}

/* 解压json格式数据 */
void coder_json::decoder(const uint8_t *in, const int64_t in_len, Json::Value &out)
{
    JSONCPP_STRING err;
    Json::CharReaderBuilder builder;
    if (in_len <= 0)
        return;

    ZSTD_DCtx *const dctx = ZSTD_createDCtx();
    check_exit(dctx != NULL, ERR_INTERNEL, "ZSTD_createDCtx() failed!");

    size_t lastRet = 0;
    std::string outdata = "";
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    uint8_t *buffOut;
    safe_alloc(buffOutSize, uint8_t, buffOut);
    ZSTD_inBuffer input = {in, (size_t)in_len, 0};
    while (input.pos < input.size)
    {
        ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
        size_t const ret = ZSTD_decompressStream(dctx, &output, &input);

        if (ZSTD_isError(ret))
        {
            free(buffOut);
            fprintf(stderr, "ZSTD_decompress error : %s \n",
                    ZSTD_getErrorName(ret));
            manage::instance().exit(ERR_INTERNEL);
        }

        outdata += std::string((char *)buffOut, output.pos);
        lastRet = ret;
    }

    if (lastRet != 0)
    {
        /* The last return value from ZSTD_decompressStream did not end on a
         * frame, but we reached the end of the file! We assume this is an
         * error, and the input was truncated.
         */
        ZSTD_freeDCtx(dctx);
        free(buffOut);
        fprintf(stderr, "EOF before end of stream: %zu\n", lastRet);
        manage::instance().exit(ERR_INTERNEL);
    }

    const std::unique_ptr<Json::CharReader> reader(builder.newCharReader());
    if (!reader->parse((char *)(outdata.c_str()), (char *)(outdata.c_str()) + outdata.size(), &out,
                       &err))
    {
        ZSTD_freeDCtx(dctx);
        free(buffOut);
        check_exit(false, ERR_INTERNEL, "json decoder failed");
    }

    ZSTD_freeDCtx(dctx);
    free(buffOut);
}