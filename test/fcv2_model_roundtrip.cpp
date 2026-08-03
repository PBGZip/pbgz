/*
 * fcv2 模型快照独立测试。
 *
 * 训练数据、码流和损坏快照均在内存中生成，避免测试依赖大数据文件；每个断言对应模型
 * 快照接口的一项公开保证，失败立即以非零状态退出，适合直接作为构建后的冒烟测试。
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <vector>

#include "coder/coder_fcv2.h"

namespace {

const size_t kHeaderBytes = 43;

size_t expectedBlobSize(int alpha)
{
    size_t prevStates = (size_t)alpha + 1;
    size_t counters = 64 + 96 * 64 + prevStates * 96 * 64 +
        prevStates * prevStates * 16 * 2 * 64;
    return kHeaderBytes + (size_t)alpha * 3 + counters * 2 + 16 * 64 * 4 * 4;
}

std::vector<uint8_t> makeQualityData()
{
    std::vector<uint8_t> data;
    data.reserve(128 * 151);
    for (int record = 0; record < 128; record++) {
        for (int cycle = 0; cycle < 151; cycle++) {
            int value = (record * 11 + cycle * 5 + cycle / 9) % 38;
            data.push_back((uint8_t)('!' + value));
        }
    }
    return data;
}

std::vector<uint32_t> frequencies(const std::vector<uint8_t>& data)
{
    std::vector<uint32_t> result(256, 0);
    for (size_t i = 0; i < data.size(); i++) result[data[i]]++;
    return result;
}

bool trainAndExport(const std::vector<uint8_t>& data, const std::vector<uint32_t>& freq,
                    std::vector<uint8_t>& blob)
{
    std::vector<uint8_t> compressed(data.size() * 2 + 4096, 0);
    coder_io io(compressed.data(), (int32_t)compressed.size());
    coder_fcv2 coder(&io, freq);
    for (size_t offset = 0; offset < data.size(); offset += 151) {
        coder.encode_record(data.data() + offset, 151, ((offset / 151) & 1) != 0);
    }
    return coder.encode_flush() > 0 && coder.export_model(blob);
}

bool rejectsAndFallsBack(const std::vector<uint8_t>& malformed,
                         const std::vector<uint32_t>& freq)
{
    std::vector<uint8_t> storage(4096, 0);
    coder_io io(storage.data(), (int32_t)storage.size());
    bool loaded = true;
    coder_fcv2 coder(&io, freq, malformed, &loaded);
    std::vector<uint8_t> fallback;
    return !loaded && coder.export_model(fallback) && !fallback.empty();
}

bool check(const char* name, bool result)
{
    printf("%s: %s\n", name, result ? "PASS" : "FAIL");
    return result;
}

} /* namespace */

int main()
{
    const std::vector<uint8_t> input = makeQualityData();
    const std::vector<uint32_t> freq = frequencies(input);
    std::vector<uint8_t> blob;

    const bool exported = trainAndExport(input, freq, blob);
    const bool sizeCorrect = exported && blob.size() == expectedBlobSize(38);
    bool allPassed = check("6a 导出非空且长度精确", sizeCorrect);

    std::vector<uint8_t> reloadStorage(4096, 0);
    coder_io reloadIo(reloadStorage.data(), (int32_t)reloadStorage.size());
    bool loaded = false;
    coder_fcv2 reloaded(&reloadIo, freq, blob, &loaded);
    std::vector<uint8_t> reexported;
    allPassed = check("6b 载入后逐字节重导出", loaded && reloaded.export_model(reexported) &&
                      reexported == blob) && allPassed;

    std::vector<uint8_t> compressed(input.size() * 2 + 4096, 0);
    coder_io encoderIo(compressed.data(), (int32_t)compressed.size());
    bool encoderLoaded = false;
    coder_fcv2 encoder(&encoderIo, freq, blob, &encoderLoaded);
    for (size_t offset = 0; offset < input.size(); offset += 151) {
        encoder.encode_record(input.data() + offset, 151, ((offset / 151) & 1) != 0);
    }
    int32_t compressedLength = encoder.encode_flush();

    coder_io decoderIo(compressed.data(), compressedLength > 0 ? compressedLength : 0);
    bool decoderLoaded = false;
    coder_fcv2 decoder(&decoderIo, freq, blob, &decoderLoaded);
    std::vector<uint8_t> decoded(input.size(), 0);
    int32_t decodeStart = decoder.begin_decode();
    bool roundTrip = encoderLoaded && decoderLoaded && compressedLength > 0 && decodeStart == 0;
    for (size_t offset = 0; roundTrip && offset < input.size(); offset += 151) {
        roundTrip = decoder.decode_record(decoded.data() + offset, 151) == 151;
    }
    allPassed = check("6c 先验加载后的编解码往返", roundTrip && decoded == input) && allPassed;

    std::vector<uint8_t> truncated(blob.begin(), blob.empty() ? blob.end() : blob.end() - 1);
    std::vector<uint8_t> badMagic = blob;
    std::vector<uint8_t> badVersion = blob;
    std::vector<uint8_t> badLength = blob;
    if (!badMagic.empty()) badMagic[0] ^= 0xFF;
    if (badVersion.size() > 10) badVersion[9]++;
    if (badLength.size() > 24) badLength[21]++;
    const bool corruptionsRejected = rejectsAndFallsBack(truncated, freq) &&
        rejectsAndFallsBack(badMagic, freq) && rejectsAndFallsBack(badVersion, freq) &&
        rejectsAndFallsBack(badLength, freq);
    allPassed = check("6d 截断、魔数、版本和长度损坏均回退", corruptionsRejected) && allPassed;

    printf("模型快照长度: alpha=38 为 %zu 字节，alpha=64 为 %zu 字节\n",
           expectedBlobSize(38), expectedBlobSize(64));
    return allPassed ? 0 : 1;
}
