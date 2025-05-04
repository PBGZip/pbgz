#ifndef _ERR_CODE_H_
#define _ERR_CODE_H_

enum errcode
{
    /* 错误的参数 */
    ERR_BAD_ARGS = (-1),

    /* 内存不足*/
    ERR_MEM_NOENOUGH = ((ERR_BAD_ARGS)-1),

    /* 内部错误 */
    ERR_INTERNEL = ((ERR_MEM_NOENOUGH)-1),

    /* 文件不存在*/
    ERR_FILE_NOTFOUND = ((ERR_INTERNEL)-1),

    /* 文件读错误 */
    ERR_FILE_READ = ((ERR_FILE_NOTFOUND)-1),

    /* 文件写错误 */
    ERR_FILE_WRITE = ((ERR_FILE_READ)-1),

    /* 文件格式错误 */
    ERR_FILE_FORMAT = ((ERR_FILE_WRITE)-1),

    /* 版本不兼容 */
    ERR_VERSION_MATCH = ((ERR_FILE_FORMAT)-1),

    /* 正常无错误 */
    ERR_NO = 0
};

#endif