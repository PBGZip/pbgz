

## 概述
PBGZ设计为一种开放、灵活的，面向超大型数据压缩的工具。主要支持一下功能：
- 支持极高的压缩率，达到业界领先水平。
- 较快的压缩和解压速度，较小的资源占用。
- 支持超大文件随机读，如指定位置的文件内容的解压。
- 支持超大文件增量更新。
- 支持大文件直接拼接，直接使用操作系统的cat命令，而非专用工具。
- 支持压缩文件在部分损毁的情况下可忽略损毁部分，继续解压。
- 支持流输入的压缩，解压和读取。比如数据来自标准输入（stdin）,不清楚数据的总大小。

	
##	PBGZ文件格式设计
PBGZ文件分为三个部分，也就是三种块：文件头块、文件描述信息块和数据块，文件头，文件描述信息在压缩一个文件时，只存在一份，数据块可以有多份。当两个压缩文件拼接时，则后一个文件的文件头和文件描述信息追加在前一个文件数据块之后。
### 文件头块包括如下信息：
- magic：固定为PBGZ，4个字节,用于标记此文件给pbgz格式的文件。
- version:分别存储版本号的major（主版本）,minor（次版本）,patch（补丁版本），每个占一个字节，共3个字节。
### 文件描述信息块包括如下信息：
- magic：0x000000FB, 4个字节，用于标记数据为文件描述信息
- origin format: 2个字节，用于记录原始文件格式，不同文件对于meta信息记录的要求不同。
- file meta length: 文件meta信息长度，4个字节。
- file meta json: json格式，文件的描述信息，字符串以’\0’结尾。
- file meta checksum：描述信息的校验和，8个字节。
### 数据块包括如下内容：
- magic：固定0x000000DB,4个字节，用于标记数据为文件数据块
- block data meta length:块数据meta信息长度，4个字节。
- block data meta: 数据块描述信息。Json格式，以’\0’结尾。
- data size: 数据块大小，4个字节
- data content：数据块内容
- data block checksum: 数据块的校验和，8个字节
- origin data checksum:原始数据的校验和，8个字节
			
### 整体文件结构图如下：
<table bold>
<tr align="center">
    <td colspan = 7> 文件头 </td>
</tr>
<tr align="center">
    <td  colspan = 3>magic</td>
    <td  colspan = 4>version</td>
</tr>
<tr align="center">
    <td colspan = 7>文件描述信息</td>
</tr>
<tr align="center">
    <td>magic</td>
    <td>origin format</td>
    <td colspan = 2>meta json length</td>
    <td>meta json</td>
    <td colspan = 2>meta json checksum</td>
</tr>
<tr align="center">
    <td colspan = 7> 文件内容块1 </td>
<tr>
<tr align="center">
    <td>magic</td>
    <td>meta data length</td>
    <td>meta data</td>
    <td>data size</td>
    <td>data content</td>
    <td>data block checksum</td>
    <td>origin data checksum</td>
</tr>
<tr align="center">
    <td colspan = 7> 文件内容块2 </td>
<tr>
<tr align="center">
    <td>magic</td>
    <td>meta data length</td>
    <td>meta data</td>
    <td>data size</td>
    <td>data content</td>
    <td>data block checksum</td>
    <td>origin data checksum</td>
</tr>
</table>

