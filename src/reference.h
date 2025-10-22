#pragma once

#include <stdint.h>
#include <string>
#include <json/json.h>
#include "utils/guard_bar.h"

typedef std::vector<std::pair<std::pair<uint32_t*, uint32_t>, uint32_t>> HashTable;
typedef std::pair<std::pair<int64_t, uint8_t*>, std::pair<int64_t, int64_t>> RefeInfo;

class Reference {
private:
    typedef struct {
        uint32_t baseGroupPos;
        uint32_t hashBucket;
    } BaseGroupHash;

public:
    Reference(const std::string& fasqaName, uint32_t threadNum);

    virtual ~Reference();

    /// @brief 创建索引表
    /// @return  bool true成功 / false 失败
    bool makeIndex();

    /// @brief 根据hash值查询参考基因组中的位置信息
    /// @param hash 输入的hash值
    /// @param length 输出参数，返回该hash bucket中位置信息的数量
    /// @return 指向位置信息数组的指针，数组包含length个uint32_t值
    const uint32_t* queryPosition(const uint32_t &hash, uint32_t &length);

    /// @brief 获取建索引时使用的碱基组长度
    /// @return 碱基组长度，固定为31
    uint32_t getBaseGroupLength() const;

    /// @brief 获取建索引时使用的碱基组步长
    /// @return 碱基组步长，固定为32
    uint32_t getBaseGroupStep()  const;

    /// @brief 获取参考基因组压缩后的数据缓冲区
    /// @return 指向压缩数据缓冲区的指针，每个字节包含4个碱基的2位编码
    const uint8_t* getSquash() const;

    /// @brief 获取参考基因组压缩数据的长度
    /// @return 压缩数据的字节长度
    int64_t getSquashLength() const;

    /// @brief 获取FASTA参考基因组文件名
    /// @return FASTA文件名的字符串引用
    const std::string& getFastaFileName() const;

    /// @brief 获取FASTA参考基因组文件的内容长度
    /// @return FASTA文件的字节长度
    int64_t getFastaLength() const;

    /// @brief 获取FASTA参考基因组文件的MD5校验和
    /// @return MD5校验和字符串
    const std::string& getFastaChecksum() const;

    /// @brief 获取NI索引文件的完整路径
    /// @return NI文件路径的字符串引用
    const std::string& getNiFilePath() const;

    /// @brief 清理参考基因组中未匹配的区域
    /// 将指定范围内未匹配的squash字节设置为0
    /// @param startPos 起始的squash位置（字节索引）
    /// @param length 需要清理的长度（字节数）
    void sanitizeRefSquash(int64_t startPos, int64_t length);

    /// @brief 更新已匹配基因区域的信息
    /// 标记指定位置范围的参考基因组区域为已匹配
    /// @param actgPos 起始的ACTG碱基位置
    /// @param matchLength 匹配的碱基长度
    void updateMatchedGene(uint64_t actgPos, uint32_t matchLength);

    /// @brief 根据参考基因组文件生成NI索引文件路径
    /// @param niFile 输出参数，返回生成的NI文件完整路径
    void getNiFileFromReference(std::string& niFile);

    /// @brief 从NI文件初始化参考基因组压缩数据
    /// 读取NI文件并初始化refGeneSquash缓冲区
    /// @return 成功返回true，失败返回false
    bool initSquashByNiFile();

    /// @brief 通过流式方式初始化参考基因组压缩数据
    /// 分配指定大小的缓冲区用于存储压缩数据
    /// @param squashLength 压缩数据的期望长度
    /// @return 指向分配的压缩数据缓冲区的指针
    uint8_t* initSquashByStream(int64_t squashLength);

    /// @brief 获取指定位置范围的ACTG碱基序列
    /// 从参考基因组中提取指定位置和长度的ACTG碱基序列
    /// @param out 输出缓冲区，用于存储ACTG字符
    /// @param outLength 输出缓冲区的长度
    /// @param actgPos 起始的ACTG碱基位置
    void getStretchActg(uint8_t* out, uint32_t outLength, uint64_t actgPos);

    /// @brief 获取指定位置范围的2位编码序列
    /// 从参考基因组中提取指定位置和长度的2位编码序列（每个字符存储2位）
    /// @param out 输出缓冲区，用于存储2位编码
    /// @param outLength 输出缓冲区的长度
    /// @param actgPos 起始的ACTG碱基位置
    void getStretch2Bits1Char(uint8_t* out, uint32_t outLength, uint64_t actgPos);

    /// @brief 将2位编码转换为ACTG字符序列
    /// 将输入的2位编码数据转换为对应的ACTG字符
    /// @param src2Bit 输入的2位编码数据
    /// @param src2BitsLen 输入数据的长度
    /// @param dstActg 输出缓冲区，用于存储ACTG字符
    void getActgFrom2Bits(const uint8_t* src2Bit, uint32_t src2BitsLen, uint8_t* dstActg);

private:
    /// @brief 检查参考基因组参数的有效性
    /// 验证baseGroupStep和baseGroupLen等关键参数是否符合要求
    /// @return 成功返回PBGZ_ERR_OK，失败返回错误码
    int32_t referencCheck();

    /// @brief 从参考基因组中提取碱基组并计算hash值
    /// 多线程处理，将参考基因组按baseGroupStep步长切分为baseGroupLen长度的碱基组，
    /// 计算每个碱基组的hash值并存储到bgHash数组中
    /// @param bgHash 输出参数，指向存储碱基组hash信息的数组
    void makeIndexFetchBaseGroup(BaseGroupHash *&bgHash);

    /// @brief 计算hash表的大小和每个bucket的元素数量
    /// 统计每个hash bucket中的元素数量，为hash表分配内存空间
    /// @param hashTable 输出参数，包含hash表大小信息的向量
    void makeIndexCalcHashTableSize(HashTable &hashTable);

    /// @brief 初始化hash表结构
    /// 根据计算得到的hash表大小信息，初始化hash表的内存布局，
    /// 设置每个bucket的起始位置和冲突处理机制
    /// @param hashTable hash表大小信息
    /// @param hashBucketCurPos 输出参数，每个bucket的当前位置指针数组
    void makeIndexInitHashTable(const HashTable& hashTable, uint32_t* &hashBucketCurPos);

    /// @brief 构建hash表内容
    /// 将碱基组的位置信息填充到hash表中，处理hash冲突
    /// @param bgHash 碱基组hash信息数组
    /// @param hashBucketCurPos 每个bucket的当前位置指针数组
    void makeIndexBuildHashTable(const BaseGroupHash* bgHash, uint32_t* &hashBucketCurPos);
    
    /// @brief 对hash表进行排序
    /// 对每个bucket中的位置信息进行排序，优化查询性能
    void makeIndexSortHashTable();

    /// @brief 检查NI索引文件是否有效
    /// 验证NI文件的存在性、格式正确性以及与参考基因组文件的匹配性
    /// @param niFile NI文件路径
    /// @return 有效返回true，无效返回false
    bool isNiFileValid(const std::string& niFile);

    /// @brief 创建NI索引文件
    /// 读取FASTA参考基因组文件，进行压缩编码，生成包含压缩数据和元数据的NI索引文件
    /// @param niFile 要创建的NI文件路径
    /// @return 成功返回true，失败返回false
    bool makeNiFile(const std::string& niFile);

    /// @brief 调试函数：导出hash表内容到文件
    /// 将hash表的内容写入"hash_table"文件，用于调试和分析
    void dumpHashTable();
    
private:
    // 进度条
    GuardBar* guardBar;
    int64_t gbCurrent;
    int64_t gbTotal;
    
    // fasta文件校验和, 当前用md5
    std::string fastaChecksum;
    // ni文件的路径
    std::string niFilePath;
    // fasqa文件内容长度
    int64_t fastaLength;
    // 建立索引表时取碱基作为key的长度，必须为奇数
    const uint32_t baseGroupLen = 31;
    // 建立索引表时碱基的步长
    const uint32_t baseGroupStep = 32;
    // 并发度
    uint32_t parallel;

    // hash table信息
    // 每个hash buckets元素个数
    uint32_t* hashBucketCnt;

    //  记录参考基因组碱基位置信息的hash table buffer
    uint32_t* hashTableBuffer;

    // 参考基因组文件
    std::string refGeneFile;

    // 参考基因actg编码后的buffer
    uint8_t* refGeneSquash;

    // 参考基因actg编码后的buffer 对应的长度
    int64_t refGeneSquashlen;

    // 参考基因组actg编码后,匹配上了base的buffer
    uint8_t* refGeneSquashMatched;

    // 对应的长度
    uint64_t refGeneSquashMatchedlen;

    // 缓存文件
    Json::Value ref2niCache;

    // actg hash信息
    const int32_t hashBuckets = (32 << 20); // 32M

    // hash buckets对应的mask
    const int32_t hashMask = hashBuckets - 1;

    // 1 byte squashed actg stretch to 4 actg
    const uint32_t actgStretch[256] = {
        0x41414141, 0x43414141, 0x54414141, 0x47414141, 0x41434141, 0x43434141, 0x54434141, 0x47434141, 0x41544141, 0x43544141, 0x54544141, 0x47544141, 0x41474141, 0x43474141, 0x54474141, 0x47474141,
        0x41414341, 0x43414341, 0x54414341, 0x47414341, 0x41434341, 0x43434341, 0x54434341, 0x47434341, 0x41544341, 0x43544341, 0x54544341, 0x47544341, 0x41474341, 0x43474341, 0x54474341, 0x47474341,
        0x41415441, 0x43415441, 0x54415441, 0x47415441, 0x41435441, 0x43435441, 0x54435441, 0x47435441, 0x41545441, 0x43545441, 0x54545441, 0x47545441, 0x41475441, 0x43475441, 0x54475441, 0x47475441,
        0x41414741, 0x43414741, 0x54414741, 0x47414741, 0x41434741, 0x43434741, 0x54434741, 0x47434741, 0x41544741, 0x43544741, 0x54544741, 0x47544741, 0x41474741, 0x43474741, 0x54474741, 0x47474741,
        0x41414143, 0x43414143, 0x54414143, 0x47414143, 0x41434143, 0x43434143, 0x54434143, 0x47434143, 0x41544143, 0x43544143, 0x54544143, 0x47544143, 0x41474143, 0x43474143, 0x54474143, 0x47474143,
        0x41414343, 0x43414343, 0x54414343, 0x47414343, 0x41434343, 0x43434343, 0x54434343, 0x47434343, 0x41544343, 0x43544343, 0x54544343, 0x47544343, 0x41474343, 0x43474343, 0x54474343, 0x47474343,
        0x41415443, 0x43415443, 0x54415443, 0x47415443, 0x41435443, 0x43435443, 0x54435443, 0x47435443, 0x41545443, 0x43545443, 0x54545443, 0x47545443, 0x41475443, 0x43475443, 0x54475443, 0x47475443,
        0x41414743, 0x43414743, 0x54414743, 0x47414743, 0x41434743, 0x43434743, 0x54434743, 0x47434743, 0x41544743, 0x43544743, 0x54544743, 0x47544743, 0x41474743, 0x43474743, 0x54474743, 0x47474743,
        0x41414154, 0x43414154, 0x54414154, 0x47414154, 0x41434154, 0x43434154, 0x54434154, 0x47434154, 0x41544154, 0x43544154, 0x54544154, 0x47544154, 0x41474154, 0x43474154, 0x54474154, 0x47474154,
        0x41414354, 0x43414354, 0x54414354, 0x47414354, 0x41434354, 0x43434354, 0x54434354, 0x47434354, 0x41544354, 0x43544354, 0x54544354, 0x47544354, 0x41474354, 0x43474354, 0x54474354, 0x47474354,
        0x41415454, 0x43415454, 0x54415454, 0x47415454, 0x41435454, 0x43435454, 0x54435454, 0x47435454, 0x41545454, 0x43545454, 0x54545454, 0x47545454, 0x41475454, 0x43475454, 0x54475454, 0x47475454,
        0x41414754, 0x43414754, 0x54414754, 0x47414754, 0x41434754, 0x43434754, 0x54434754, 0x47434754, 0x41544754, 0x43544754, 0x54544754, 0x47544754, 0x41474754, 0x43474754, 0x54474754, 0x47474754,
        0x41414147, 0x43414147, 0x54414147, 0x47414147, 0x41434147, 0x43434147, 0x54434147, 0x47434147, 0x41544147, 0x43544147, 0x54544147, 0x47544147, 0x41474147, 0x43474147, 0x54474147, 0x47474147,
        0x41414347, 0x43414347, 0x54414347, 0x47414347, 0x41434347, 0x43434347, 0x54434347, 0x47434347, 0x41544347, 0x43544347, 0x54544347, 0x47544347, 0x41474347, 0x43474347, 0x54474347, 0x47474347,
        0x41415447, 0x43415447, 0x54415447, 0x47415447, 0x41435447, 0x43435447, 0x54435447, 0x47435447, 0x41545447, 0x43545447, 0x54545447, 0x47545447, 0x41475447, 0x43475447, 0x54475447, 0x47475447,
        0x41414747, 0x43414747, 0x54414747, 0x47414747, 0x41434747, 0x43434747, 0x54434747, 0x47434747, 0x41544747, 0x43544747, 0x54544747, 0x47544747, 0x41474747, 0x43474747, 0x54474747, 0x47474747};

    // 1 byte squashed actg stretch to 4 2bits
    const uint32_t actgStretch2bits[256] = {
        0x00000000, 0x01000000, 0x02000000, 0x03000000, 0x00010000, 0x01010000, 0x02010000, 0x03010000, 0x00020000, 0x01020000, 0x02020000, 0x03020000, 0x00030000, 0x01030000, 0x02030000, 0x03030000,
        0x00000100, 0x01000100, 0x02000100, 0x03000100, 0x00010100, 0x01010100, 0x02010100, 0x03010100, 0x00020100, 0x01020100, 0x02020100, 0x03020100, 0x00030100, 0x01030100, 0x02030100, 0x03030100,
        0x00000200, 0x01000200, 0x02000200, 0x03000200, 0x00010200, 0x01010200, 0x02010200, 0x03010200, 0x00020200, 0x01020200, 0x02020200, 0x03020200, 0x00030200, 0x01030200, 0x02030200, 0x03030200,
        0x00000300, 0x01000300, 0x02000300, 0x03000300, 0x00010300, 0x01010300, 0x02010300, 0x03010300, 0x00020300, 0x01020300, 0x02020300, 0x03020300, 0x00030300, 0x01030300, 0x02030300, 0x03030300,
        0x00000001, 0x01000001, 0x02000001, 0x03000001, 0x00010001, 0x01010001, 0x02010001, 0x03010001, 0x00020001, 0x01020001, 0x02020001, 0x03020001, 0x00030001, 0x01030001, 0x02030001, 0x03030001,
        0x00000101, 0x01000101, 0x02000101, 0x03000101, 0x00010101, 0x01010101, 0x02010101, 0x03010101, 0x00020101, 0x01020101, 0x02020101, 0x03020101, 0x00030101, 0x01030101, 0x02030101, 0x03030101,
        0x00000201, 0x01000201, 0x02000201, 0x03000201, 0x00010201, 0x01010201, 0x02010201, 0x03010201, 0x00020201, 0x01020201, 0x02020201, 0x03020201, 0x00030201, 0x01030201, 0x02030201, 0x03030201,
        0x00000301, 0x01000301, 0x02000301, 0x03000301, 0x00010301, 0x01010301, 0x02010301, 0x03010301, 0x00020301, 0x01020301, 0x02020301, 0x03020301, 0x00030301, 0x01030301, 0x02030301, 0x03030301,
        0x00000002, 0x01000002, 0x02000002, 0x03000002, 0x00010002, 0x01010002, 0x02010002, 0x03010002, 0x00020002, 0x01020002, 0x02020002, 0x03020002, 0x00030002, 0x01030002, 0x02030002, 0x03030002,
        0x00000102, 0x01000102, 0x02000102, 0x03000102, 0x00010102, 0x01010102, 0x02010102, 0x03010102, 0x00020102, 0x01020102, 0x02020102, 0x03020102, 0x00030102, 0x01030102, 0x02030102, 0x03030102,
        0x00000202, 0x01000202, 0x02000202, 0x03000202, 0x00010202, 0x01010202, 0x02010202, 0x03010202, 0x00020202, 0x01020202, 0x02020202, 0x03020202, 0x00030202, 0x01030202, 0x02030202, 0x03030202,
        0x00000302, 0x01000302, 0x02000302, 0x03000302, 0x00010302, 0x01010302, 0x02010302, 0x03010302, 0x00020302, 0x01020302, 0x02020302, 0x03020302, 0x00030302, 0x01030302, 0x02030302, 0x03030302,
        0x00000003, 0x01000003, 0x02000003, 0x03000003, 0x00010003, 0x01010003, 0x02010003, 0x03010003, 0x00020003, 0x01020003, 0x02020003, 0x03020003, 0x00030003, 0x01030003, 0x02030003, 0x03030003,
        0x00000103, 0x01000103, 0x02000103, 0x03000103, 0x00010103, 0x01010103, 0x02010103, 0x03010103, 0x00020103, 0x01020103, 0x02020103, 0x03020103, 0x00030103, 0x01030103, 0x02030103, 0x03030103,
        0x00000203, 0x01000203, 0x02000203, 0x03000203, 0x00010203, 0x01010203, 0x02010203, 0x03010203, 0x00020203, 0x01020203, 0x02020203, 0x03020203, 0x00030203, 0x01030203, 0x02030203, 0x03030203,
        0x00000303, 0x01000303, 0x02000303, 0x03000303, 0x00010303, 0x01010303, 0x02010303, 0x03010303, 0x00020303, 0x01020303, 0x02020303, 0x03020303, 0x00030303, 0x01030303, 0x02030303, 0x03030303};

};