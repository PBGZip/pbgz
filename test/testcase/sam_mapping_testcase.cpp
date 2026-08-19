/*
 * sam_mapping_testcase.cpp - Test cases for SAM mapping functionality
 * Copyright (c) 2021, Baidu.com, Inc.
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>

#define private public
#include "sam_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include <actg.h>
#undef private

namespace SamMappingData {
    const std::string testSamFile = "test.sam";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamMappingTest : public ::testing::Test {
public:
    void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);

		pInBlock = new RoughIOBlock(SamMappingData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(SamMappingData::MAX_BLOCK_SIZE);

        generateSamFile(SamMappingData::testSamFile);

        ConfigManager::getInstance().logLevel = LogLevel::WARNING;
    }

    void TearDown() override {
        if (pInBlock != nullptr) {
			delete pInBlock;
            pInBlock = nullptr;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
            pOutBlock = nullptr;
        }
        std::remove(SamMappingData::testSamFile.c_str()); 
    }

    void generateSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            // If unable to open in current directory, try to create in test directory
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                // If still fails, try to create in build directory
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }
        
        if (!file.is_open()) {
            return; // Unable to create file
        }

        file << "@SQ\tSN:CM000663.2\tLN:248956422\n"
             << "@SQ\tSN:KI270706.1\tLN:175055\n"
             << "@SQ\tSN:KI270707.1\tLN:32032\n"
             << "@SQ\tSN:KI270708.1\tLN:127682\n"
             << "@SQ\tSN:KI270709.1\tLN:66860\n"
             << "@SQ\tSN:KI270710.1\tLN:40176\n"
             << "@SQ\tSN:KI270711.1\tLN:42210\n"
             << "@SQ\tSN:KI270712.1\tLN:176043\n"
             << "@SQ\tSN:KI270713.1\tLN:40745\n"
             << "@SQ\tSN:KI270714.1\tLN:41717\n"
             << "@SQ\tSN:CM000664.2\tLN:242193529\n"
             << "@SQ\tSN:KI270715.1\tLN:161471\n"
             << "@SQ\tSN:KI270716.1\tLN:153799\n"
             << "@SQ\tSN:CM000665.2\tLN:198295559\n"
             << "@SQ\tSN:GL000221.1\tLN:155397\n"
             << "@SQ\tSN:CM000666.2\tLN:190214555\n"
             << "@SQ\tSN:GL000008.2\tLN:209709\n"
             << "@SQ\tSN:CM000667.2\tLN:181538259\n"
             << "@SQ\tSN:GL000208.1\tLN:92689\n"
             << "@SQ\tSN:CM000668.2\tLN:170805979\n"
             << "@SQ\tSN:CM000669.2\tLN:159345973\n"
             << "@SQ\tSN:CM000670.2\tLN:145138636\n"
             << "@SQ\tSN:CM000671.2\tLN:138394717\n"
             << "@SQ\tSN:KI270717.1\tLN:40062\n"
             << "@SQ\tSN:KI270718.1\tLN:38054\n"
             << "@SQ\tSN:KI270719.1\tLN:176845\n"
             << "@SQ\tSN:KI270720.1\tLN:39050\n"
             << "@SQ\tSN:CM000672.2\tLN:133797422\n"
             << "@SQ\tSN:CM000673.2\tLN:135086622\n"
             << "@SQ\tSN:KI270721.1\tLN:100316\n"
             << "@SQ\tSN:CM000674.2\tLN:133275309\n"
             << "@SQ\tSN:CM000675.2\tLN:114364328\n"
             << "@SQ\tSN:CM000676.2\tLN:107043718\n"
             << "@SQ\tSN:GL000009.2\tLN:201709\n"
             << "@SQ\tSN:GL000225.1\tLN:211173\n"
             << "@SQ\tSN:KI270722.1\tLN:194050\n"
             << "@SQ\tSN:GL000194.1\tLN:191469\n"
             << "@SQ\tSN:KI270723.1\tLN:38115\n"
             << "@SQ\tSN:KI270724.1\tLN:39555\n"
             << "@SQ\tSN:KI270725.1\tLN:172810\n"
             << "@SQ\tSN:KI270726.1\tLN:43739\n"
             << "@SQ\tSN:CM000677.2\tLN:101991189\n"
             << "@SQ\tSN:KI270727.1\tLN:448248\n"
             << "@SQ\tSN:CM000678.2\tLN:90338345\n"
             << "@SQ\tSN:KI270728.1\tLN:1872759\n"
             << "@SQ\tSN:CM000679.2\tLN:83257441\n"
             << "@SQ\tSN:GL000205.2\tLN:185591\n"
             << "@SQ\tSN:KI270729.1\tLN:280839\n"
             << "@SQ\tSN:KI270730.1\tLN:112551\n"
             << "@SQ\tSN:CM000680.2\tLN:80373285\n"
             << "@SQ\tSN:CM000681.2\tLN:58617616\n"
             << "@SQ\tSN:CM000682.2\tLN:64444167\n"
             << "@SQ\tSN:CM000683.2\tLN:46709983\n"
             << "@SQ\tSN:CM000684.2\tLN:50818468\n"
             << "@SQ\tSN:KI270731.1\tLN:150754\n"
             << "@SQ\tSN:KI270732.1\tLN:41543\n"
             << "@SQ\tSN:KI270733.1\tLN:179772\n"
             << "@SQ\tSN:KI270734.1\tLN:165050\n"
             << "@SQ\tSN:KI270735.1\tLN:42811\n"
             << "@SQ\tSN:KI270736.1\tLN:181920\n"
             << "@SQ\tSN:KI270737.1\tLN:103838\n"
             << "@SQ\tSN:KI270738.1\tLN:99375\n"
             << "@SQ\tSN:KI270739.1\tLN:73985\n"
             << "@SQ\tSN:CM000685.2\tLN:156040895\n"
             << "@SQ\tSN:CM000686.2\tLN:57227415\n"
             << "@SQ\tSN:KI270740.1\tLN:37240\n"
             << "@SQ\tSN:KI270302.1\tLN:2274\n"
             << "@SQ\tSN:KI270304.1\tLN:2165\n"
             << "@SQ\tSN:KI270303.1\tLN:1942\n"
             << "@SQ\tSN:KI270305.1\tLN:1472\n"
             << "@SQ\tSN:KI270322.1\tLN:21476\n"
             << "@SQ\tSN:KI270320.1\tLN:4416\n"
             << "@SQ\tSN:KI270310.1\tLN:1201\n"
             << "@SQ\tSN:KI270316.1\tLN:1444\n"
             << "@SQ\tSN:KI270315.1\tLN:2276\n"
             << "@SQ\tSN:KI270312.1\tLN:998\n"
             << "@SQ\tSN:KI270311.1\tLN:12399\n"
             << "@SQ\tSN:KI270317.1\tLN:37690\n"
             << "@SQ\tSN:KI270412.1\tLN:1179\n"
             << "@SQ\tSN:KI270411.1\tLN:2646\n"
             << "@SQ\tSN:KI270414.1\tLN:2489\n"
             << "@SQ\tSN:KI270419.1\tLN:1029\n"
             << "@SQ\tSN:KI270418.1\tLN:2145\n"
             << "@SQ\tSN:KI270420.1\tLN:2321\n"
             << "@SQ\tSN:KI270424.1\tLN:2140\n"
             << "@SQ\tSN:KI270417.1\tLN:2043\n"
             << "@SQ\tSN:KI270422.1\tLN:1445\n"
             << "@SQ\tSN:KI270423.1\tLN:981\n"
             << "@SQ\tSN:KI270425.1\tLN:1884\n"
             << "@SQ\tSN:KI270429.1\tLN:1361\n"
             << "@SQ\tSN:KI270442.1\tLN:392061\n"
             << "@SQ\tSN:KI270466.1\tLN:1233\n"
             << "@SQ\tSN:KI270465.1\tLN:1774\n"
             << "@SQ\tSN:KI270467.1\tLN:3920\n"
             << "@SQ\tSN:KI270435.1\tLN:92983\n"
             << "@SQ\tSN:KI270438.1\tLN:112505\n"
             << "@SQ\tSN:KI270468.1\tLN:4055\n"
             << "@SQ\tSN:KI270510.1\tLN:2415\n"
             << "@SQ\tSN:KI270509.1\tLN:2318\n"
             << "@SQ\tSN:KI270518.1\tLN:2186\n"
             << "@SQ\tSN:KI270508.1\tLN:1951\n"
             << "@SQ\tSN:KI270516.1\tLN:1300\n"
             << "@SQ\tSN:KI270512.1\tLN:22689\n"
             << "@SQ\tSN:KI270519.1\tLN:138126\n"
             << "@SQ\tSN:KI270522.1\tLN:5674\n"
             << "@SQ\tSN:KI270511.1\tLN:8127\n"
             << "@SQ\tSN:KI270515.1\tLN:6361\n"
             << "@SQ\tSN:KI270507.1\tLN:5353\n"
             << "@SQ\tSN:KI270517.1\tLN:3253\n"
             << "@SQ\tSN:KI270529.1\tLN:1899\n"
             << "@SQ\tSN:KI270528.1\tLN:2983\n"
             << "@SQ\tSN:KI270530.1\tLN:2168\n"
             << "@SQ\tSN:KI270539.1\tLN:993\n"
             << "@SQ\tSN:KI270538.1\tLN:91309\n"
             << "@SQ\tSN:KI270544.1\tLN:1202\n"
             << "@SQ\tSN:KI270548.1\tLN:1599\n"
             << "@SQ\tSN:KI270583.1\tLN:1400\n"
             << "@SQ\tSN:KI270587.1\tLN:2969\n"
             << "@SQ\tSN:KI270580.1\tLN:1553\n"
             << "@SQ\tSN:KI270581.1\tLN:7046\n"
             << "@SQ\tSN:KI270579.1\tLN:31033\n"
             << "@SQ\tSN:KI270589.1\tLN:44474\n"
             << "@SQ\tSN:KI270590.1\tLN:4685\n"
             << "@SQ\tSN:KI270584.1\tLN:4513\n"
             << "@SQ\tSN:KI270582.1\tLN:6504\n"
             << "@SQ\tSN:KI270588.1\tLN:6158\n"
             << "@SQ\tSN:KI270593.1\tLN:3041\n"
             << "@SQ\tSN:KI270591.1\tLN:5796\n"
             << "@SQ\tSN:KI270330.1\tLN:1652\n"
             << "@SQ\tSN:KI270329.1\tLN:1040\n"
             << "@SQ\tSN:KI270334.1\tLN:1368\n"
             << "@SQ\tSN:KI270333.1\tLN:2699\n"
             << "@SQ\tSN:KI270335.1\tLN:1048\n"
             << "@SQ\tSN:KI270338.1\tLN:1428\n"
             << "@SQ\tSN:KI270340.1\tLN:1428\n"
             << "@SQ\tSN:KI270336.1\tLN:1026\n"
             << "@SQ\tSN:KI270337.1\tLN:1121\n"
             << "@SQ\tSN:KI270363.1\tLN:1803\n"
             << "@SQ\tSN:KI270364.1\tLN:2855\n"
             << "@SQ\tSN:KI270362.1\tLN:3530\n"
             << "@SQ\tSN:KI270366.1\tLN:8320\n"
             << "@SQ\tSN:KI270378.1\tLN:1048\n"
             << "@SQ\tSN:KI270379.1\tLN:1045\n"
             << "@SQ\tSN:KI270389.1\tLN:1298\n"
             << "@SQ\tSN:KI270390.1\tLN:2387\n"
             << "@SQ\tSN:KI270387.1\tLN:1537\n"
             << "@SQ\tSN:KI270395.1\tLN:1143\n"
             << "@SQ\tSN:KI270396.1\tLN:1880\n"
             << "@SQ\tSN:KI270388.1\tLN:1216\n"
             << "@SQ\tSN:KI270394.1\tLN:970\n"
             << "@SQ\tSN:KI270386.1\tLN:1788\n"
             << "@SQ\tSN:KI270391.1\tLN:1484\n"
             << "@SQ\tSN:KI270383.1\tLN:1750\n"
             << "@SQ\tSN:KI270393.1\tLN:1308\n"
             << "@SQ\tSN:KI270384.1\tLN:1658\n"
             << "@SQ\tSN:KI270392.1\tLN:971\n"
             << "@SQ\tSN:KI270381.1\tLN:1930\n"
             << "@SQ\tSN:KI270385.1\tLN:990\n"
             << "@SQ\tSN:KI270382.1\tLN:4215\n"
             << "@SQ\tSN:KI270376.1\tLN:1136\n"
             << "@SQ\tSN:KI270374.1\tLN:2656\n"
             << "@SQ\tSN:KI270372.1\tLN:1650\n"
             << "@SQ\tSN:KI270373.1\tLN:1451\n"
             << "@SQ\tSN:KI270375.1\tLN:2378\n"
             << "@SQ\tSN:KI270371.1\tLN:2805\n"
             << "@SQ\tSN:KI270448.1\tLN:7992\n"
             << "@SQ\tSN:KI270521.1\tLN:7642\n"
             << "@SQ\tSN:GL000195.1\tLN:182896\n"
             << "@SQ\tSN:GL000219.1\tLN:179198\n"
             << "@SQ\tSN:GL000220.1\tLN:161802\n"
             << "@SQ\tSN:GL000224.1\tLN:179693\n"
             << "@SQ\tSN:KI270741.1\tLN:157432\n"
             << "@SQ\tSN:GL000226.1\tLN:15008\n"
             << "@SQ\tSN:GL000213.1\tLN:164239\n"
             << "@SQ\tSN:KI270743.1\tLN:210658\n"
             << "@SQ\tSN:KI270744.1\tLN:168472\n"
             << "@SQ\tSN:KI270745.1\tLN:41891\n"
             << "@SQ\tSN:KI270746.1\tLN:66486\n"
             << "@SQ\tSN:KI270747.1\tLN:198735\n"
             << "@SQ\tSN:KI270748.1\tLN:93321\n"
             << "@SQ\tSN:KI270749.1\tLN:158759\n"
             << "@SQ\tSN:KI270750.1\tLN:148850\n"
             << "@SQ\tSN:KI270751.1\tLN:150742\n"
             << "@SQ\tSN:KI270752.1\tLN:27745\n"
             << "@SQ\tSN:KI270753.1\tLN:62944\n"
             << "@SQ\tSN:KI270754.1\tLN:40191\n"
             << "@SQ\tSN:KI270755.1\tLN:36723\n"
             << "@SQ\tSN:KI270756.1\tLN:79590\n"
             << "@SQ\tSN:KI270757.1\tLN:71251\n"
             << "@SQ\tSN:GL000214.1\tLN:137718\n"
             << "@SQ\tSN:KI270742.1\tLN:186739\n"
             << "@SQ\tSN:GL000216.2\tLN:176608\n"
             << "@SQ\tSN:GL000218.1\tLN:161147\n"
             << "@SQ\tSN:KQ031383.1\tLN:467143\n"
             << "@SQ\tSN:KQ983255.1\tLN:278659\n"
             << "@SQ\tSN:MU273333.1\tLN:1572686\n"
             << "@SQ\tSN:KN538361.1\tLN:305542\n"
             << "@SQ\tSN:KQ458383.1\tLN:349938\n"
             << "@SQ\tSN:KN196473.1\tLN:166200\n"
             << "@SQ\tSN:KZ208904.1\tLN:166136\n"
             << "@SQ\tSN:KZ559100.1\tLN:44955\n"
             << "@SQ\tSN:KN196472.1\tLN:186494\n"
             << "@SQ\tSN:KZ208905.1\tLN:140355\n"
             << "@SQ\tSN:KQ458382.1\tLN:141019\n"
             << "@SQ\tSN:KV880763.1\tLN:551020\n"
             << "@SQ\tSN:KN196474.1\tLN:122022\n"
             << "@SQ\tSN:MU273330.1\tLN:516764\n"
             << "@SQ\tSN:MU273335.1\tLN:211934\n"
             << "@SQ\tSN:MU273336.1\tLN:250447\n"
             << "@SQ\tSN:MU273331.1\tLN:847441\n"
             << "@SQ\tSN:KN538360.1\tLN:460100\n"
             << "@SQ\tSN:KZ208906.1\tLN:330031\n"
             << "@SQ\tSN:KQ458384.1\tLN:212205\n"
             << "@SQ\tSN:MU273334.1\tLN:210426\n"
             << "@SQ\tSN:MU273332.1\tLN:335159\n"
             << "@SQ\tSN:ML143342.1\tLN:84043\n"
             << "@SQ\tSN:MU273344.1\tLN:244725\n"
             << "@SQ\tSN:MU273345.1\tLN:174385\n"
             << "@SQ\tSN:MU273343.1\tLN:489404\n"
             << "@SQ\tSN:MU273340.1\tLN:284971\n"
             << "@SQ\tSN:KQ031384.1\tLN:481245\n"
             << "@SQ\tSN:MU273337.1\tLN:431782\n"
             << "@SQ\tSN:MU273342.1\tLN:955087\n"
             << "@SQ\tSN:KZ208907.1\tLN:181658\n"
             << "@SQ\tSN:MU273339.1\tLN:500581\n"
             << "@SQ\tSN:MU273338.1\tLN:535251\n"
             << "@SQ\tSN:MU273341.1\tLN:120381\n"
             << "@SQ\tSN:KQ983256.1\tLN:535088\n"
             << "@SQ\tSN:KZ208908.1\tLN:140361\n"
             << "@SQ\tSN:KN538363.1\tLN:365499\n"
             << "@SQ\tSN:ML143341.1\tLN:145975\n"
             << "@SQ\tSN:KN538362.1\tLN:208149\n"
             << "@SQ\tSN:KV766192.1\tLN:411654\n"
             << "@SQ\tSN:MU273347.1\tLN:301310\n"
             << "@SQ\tSN:MU273348.1\tLN:475876\n"
             << "@SQ\tSN:KN196475.1\tLN:451168\n"
             << "@SQ\tSN:KQ031385.1\tLN:373699\n"
             << "@SQ\tSN:KN538364.1\tLN:415308\n"
             << "@SQ\tSN:ML143343.1\tLN:215443\n"
             << "@SQ\tSN:KZ208909.1\tLN:175849\n"
             << "@SQ\tSN:KQ031386.1\tLN:165718\n"
             << "@SQ\tSN:KN196476.1\tLN:305979\n"
             << "@SQ\tSN:KZ559104.1\tLN:105527\n"
             << "@SQ\tSN:KZ559105.1\tLN:195063\n"
             << "@SQ\tSN:KZ559103.1\tLN:302885\n"
             << "@SQ\tSN:KZ559102.1\tLN:197752\n"
             << "@SQ\tSN:KZ559101.1\tLN:164041\n"
             << "@SQ\tSN:MU273346.1\tLN:469342\n"
             << "@SQ\tSN:ML143344.1\tLN:235734\n"
             << "@SQ\tSN:ML143347.1\tLN:176674\n"
             << "@SQ\tSN:KQ090013.1\tLN:90922\n"
             << "@SQ\tSN:MU273351.1\tLN:205691\n"
             << "@SQ\tSN:ML143346.1\tLN:53476\n"
             << "@SQ\tSN:ML143348.1\tLN:125549\n"
             << "@SQ\tSN:ML143345.1\tLN:341066\n"
             << "@SQ\tSN:ML143349.1\tLN:276109\n"
             << "@SQ\tSN:KQ090014.1\tLN:163749\n"
             << "@SQ\tSN:KQ090015.1\tLN:236512\n"
             << "@SQ\tSN:MU273349.1\tLN:308682\n"
             << "@SQ\tSN:KV766193.1\tLN:420675\n"
             << "@SQ\tSN:MU273350.1\tLN:113364\n"
             << "@SQ\tSN:KQ983258.1\tLN:205407\n"
             << "@SQ\tSN:KQ983257.1\tLN:230434\n"
             << "@SQ\tSN:MU273354.1\tLN:2101585\n"
             << "@SQ\tSN:KZ208910.1\tLN:135987\n"
             << "@SQ\tSN:MU273353.1\tLN:208405\n"
             << "@SQ\tSN:KN196477.1\tLN:139087\n"
             << "@SQ\tSN:KV575243.1\tLN:362221\n"
             << "@SQ\tSN:MU273356.1\tLN:302485\n"
             << "@SQ\tSN:ML143350.1\tLN:89956\n"
             << "@SQ\tSN:MU273355.1\tLN:508332\n"
             << "@SQ\tSN:KV575244.1\tLN:673059\n"
             << "@SQ\tSN:MU273352.1\tLN:34400\n"
             << "@SQ\tSN:KZ208911.1\tLN:242796\n"
             << "@SQ\tSN:KQ090017.1\tLN:82315\n"
             << "@SQ\tSN:MU273357.1\tLN:383128\n"
             << "@SQ\tSN:KQ031387.1\tLN:320750\n"
             << "@SQ\tSN:KN196478.1\tLN:268330\n"
             << "@SQ\tSN:KQ090016.1\tLN:245716\n"
             << "@SQ\tSN:KV766194.1\tLN:139427\n"
             << "@SQ\tSN:ML143351.1\tLN:73265\n"
             << "@SQ\tSN:ML143352.1\tLN:254759\n"
             << "@SQ\tSN:KZ559106.1\tLN:172555\n"
             << "@SQ\tSN:MU273358.1\tLN:464417\n"
             << "@SQ\tSN:KV880764.1\tLN:142129\n"
             << "@SQ\tSN:KV880765.1\tLN:468267\n"
             << "@SQ\tSN:KZ208912.1\tLN:589656\n"
             << "@SQ\tSN:KZ208913.1\tLN:680662\n"
             << "@SQ\tSN:KQ031388.1\tLN:179932\n"
             << "@SQ\tSN:MU273362.1\tLN:429744\n"
             << "@SQ\tSN:KZ208915.1\tLN:6367528\n"
             << "@SQ\tSN:KV880767.1\tLN:265876\n"
             << "@SQ\tSN:MU273359.1\tLN:150302\n"
             << "@SQ\tSN:KV880766.1\tLN:156998\n"
             << "@SQ\tSN:MU273361.1\tLN:106905\n"
             << "@SQ\tSN:KZ559107.1\tLN:103072\n"
             << "@SQ\tSN:MU273363.1\tLN:207371\n"
             << "@SQ\tSN:KZ208914.1\tLN:165120\n"
             << "@SQ\tSN:MU273360.1\tLN:39290\n"
             << "@SQ\tSN:KQ090018.1\tLN:163882\n"
             << "@SQ\tSN:KQ090019.1\tLN:134099\n"
             << "@SQ\tSN:MU273366.1\tLN:569668\n"
             << "@SQ\tSN:MU273364.1\tLN:340717\n"
             << "@SQ\tSN:MU273365.1\tLN:482250\n"
             << "@SQ\tSN:KN196479.1\tLN:330164\n"
             << "@SQ\tSN:ML143353.1\tLN:25408\n"
             << "@SQ\tSN:ML143354.1\tLN:454963\n"
             << "@SQ\tSN:KN538367.1\tLN:420164\n"
             << "@SQ\tSN:ML143355.1\tLN:292944\n"
             << "@SQ\tSN:KQ090020.1\tLN:185507\n"
             << "@SQ\tSN:KN196480.1\tLN:277797\n"
             << "@SQ\tSN:KQ090021.1\tLN:264545\n"
             << "@SQ\tSN:MU273367.1\tLN:196262\n"
             << "@SQ\tSN:KN538366.1\tLN:85284\n"
             << "@SQ\tSN:KN538365.1\tLN:14347\n"
             << "@SQ\tSN:KQ759759.2\tLN:204999\n"
             << "@SQ\tSN:MU273369.1\tLN:434831\n"
             << "@SQ\tSN:ML143358.1\tLN:270122\n"
             << "@SQ\tSN:MU273371.1\tLN:122722\n"
             << "@SQ\tSN:KN538368.1\tLN:203552\n"
             << "@SQ\tSN:ML143360.1\tLN:170928\n"
             << "@SQ\tSN:KZ559109.1\tLN:279644\n"
             << "@SQ\tSN:KZ559108.1\tLN:305244\n"
             << "@SQ\tSN:KV766195.1\tLN:140877\n"
             << "@SQ\tSN:MU273370.1\tLN:344606\n"
             << "@SQ\tSN:KZ559111.1\tLN:181167\n"
             << "@SQ\tSN:ML143359.1\tLN:217075\n"
             << "@SQ\tSN:ML143357.1\tLN:165419\n"
             << "@SQ\tSN:KZ559110.1\tLN:301637\n"
             << "@SQ\tSN:KQ090022.1\tLN:181958\n"
             << "@SQ\tSN:ML143356.1\tLN:45257\n"
             << "@SQ\tSN:MU273368.1\tLN:261194\n"
             << "@SQ\tSN:KN196481.1\tLN:108875\n"
             << "@SQ\tSN:KQ090023.1\tLN:109323\n"
             << "@SQ\tSN:KZ208916.1\tLN:1046838\n"
             << "@SQ\tSN:ML143362.1\tLN:192531\n"
             << "@SQ\tSN:KN538369.1\tLN:541038\n"
             << "@SQ\tSN:KN196482.1\tLN:211377\n"
             << "@SQ\tSN:MU273372.1\tLN:104537\n"
             << "@SQ\tSN:KZ208918.1\tLN:174808\n"
             << "@SQ\tSN:KQ759760.1\tLN:315610\n"
             << "@SQ\tSN:KZ208917.1\tLN:64689\n"
             << "@SQ\tSN:KN538370.1\tLN:86533\n"
             << "@SQ\tSN:KZ559112.1\tLN:154139\n"
             << "@SQ\tSN:ML143361.1\tLN:297568\n"
             << "@SQ\tSN:ML143366.1\tLN:409912\n"
             << "@SQ\tSN:KN538372.1\tLN:356766\n"
             << "@SQ\tSN:KQ090024.1\tLN:168146\n"
             << "@SQ\tSN:ML143363.1\tLN:7309\n"
             << "@SQ\tSN:KN196483.1\tLN:35455\n"
             << "@SQ\tSN:KN538373.1\tLN:148762\n"
             << "@SQ\tSN:KQ090025.1\tLN:123480\n"
             << "@SQ\tSN:ML143364.1\tLN:158944\n"
             << "@SQ\tSN:ML143365.1\tLN:65394\n"
             << "@SQ\tSN:KN538371.1\tLN:206320\n"
             << "@SQ\tSN:ML143367.1\tLN:399183\n"
             << "@SQ\tSN:MU273373.1\tLN:722645\n"
             << "@SQ\tSN:KZ208920.1\tLN:690932\n"
             << "@SQ\tSN:ML143368.1\tLN:264228\n"
             << "@SQ\tSN:KZ208919.1\tLN:171798\n"
             << "@SQ\tSN:KN538374.1\tLN:4998962\n"
             << "@SQ\tSN:ML143372.1\tLN:396515\n"
             << "@SQ\tSN:ML143371.1\tLN:5500449\n"
             << "@SQ\tSN:KQ031389.1\tLN:2365364\n"
             << "@SQ\tSN:MU273375.1\tLN:204007\n"
             << "@SQ\tSN:ML143370.1\tLN:369264\n"
             << "@SQ\tSN:MU273374.1\tLN:1154574\n"
             << "@SQ\tSN:ML143369.1\tLN:97763\n"
             << "@SQ\tSN:MU273376.1\tLN:87715\n"
             << "@SQ\tSN:KQ090026.1\tLN:59016\n"
             << "@SQ\tSN:KZ559113.1\tLN:480415\n"
             << "@SQ\tSN:KV880768.1\tLN:1927115\n"
             << "@SQ\tSN:ML143373.1\tLN:270967\n"
             << "@SQ\tSN:KQ090027.1\tLN:267463\n"
             << "@SQ\tSN:KZ208921.1\tLN:78609\n"
             << "@SQ\tSN:KQ031390.1\tLN:169136\n"
             << "@SQ\tSN:MU273377.1\tLN:334997\n"
             << "@SQ\tSN:KV766196.1\tLN:281919\n"
             << "@SQ\tSN:ML143374.1\tLN:137908\n"
             << "@SQ\tSN:KV575245.1\tLN:154723\n"
             << "@SQ\tSN:KV766198.1\tLN:276292\n"
             << "@SQ\tSN:MU273380.1\tLN:538541\n"
             << "@SQ\tSN:MU273378.1\tLN:372839\n"
             << "@SQ\tSN:KV766197.1\tLN:246895\n"
             << "@SQ\tSN:KZ559114.1\tLN:116753\n"
             << "@SQ\tSN:MU273383.1\tLN:172609\n"
             << "@SQ\tSN:MU273379.1\tLN:234878\n"
             << "@SQ\tSN:MU273382.1\tLN:187626\n"
             << "@SQ\tSN:ML143375.1\tLN:56695\n"
             << "@SQ\tSN:MU273381.1\tLN:144689\n"
             << "@SQ\tSN:KZ559116.1\tLN:163186\n"
             << "@SQ\tSN:KQ458385.1\tLN:205101\n"
             << "@SQ\tSN:KZ559115.1\tLN:230843\n"
             << "@SQ\tSN:KQ090028.1\tLN:407387\n"
             << "@SQ\tSN:KZ208922.1\tLN:93070\n"
             << "@SQ\tSN:MU273386.1\tLN:226166\n"
             << "@SQ\tSN:MU273387.1\tLN:89211\n"
             << "@SQ\tSN:MU273384.1\tLN:333754\n"
             << "@SQ\tSN:ML143376.1\tLN:493165\n"
             << "@SQ\tSN:KQ458386.1\tLN:405389\n"
             << "@SQ\tSN:MU273385.1\tLN:137818\n"
             << "@SQ\tSN:KN196484.1\tLN:370917\n"
             << "@SQ\tSN:KV575246.1\tLN:163926\n"
             << "@SQ\tSN:KV575247.1\tLN:170206\n"
             << "@SQ\tSN:KV575248.1\tLN:168131\n"
             << "@SQ\tSN:KV575249.1\tLN:293522\n"
             << "@SQ\tSN:KV575250.1\tLN:241058\n"
             << "@SQ\tSN:KV575251.1\tLN:159285\n"
             << "@SQ\tSN:KV575252.1\tLN:178197\n"
             << "@SQ\tSN:KV575253.1\tLN:166713\n"
             << "@SQ\tSN:KV575254.1\tLN:99845\n"
             << "@SQ\tSN:KV575255.1\tLN:161095\n"
             << "@SQ\tSN:KV575256.1\tLN:223118\n"
             << "@SQ\tSN:KV575257.1\tLN:100553\n"
             << "@SQ\tSN:KV575259.1\tLN:171263\n"
             << "@SQ\tSN:KV575260.1\tLN:145691\n"
             << "@SQ\tSN:KV575258.1\tLN:156965\n"
             << "@SQ\tSN:MU273388.1\tLN:273725\n"
             << "@SQ\tSN:MU273389.1\tLN:355731\n"
             << "@SQ\tSN:ML143377.1\tLN:519485\n"
             << "@SQ\tSN:MU273390.1\tLN:336752\n"
             << "@SQ\tSN:MU273391.1\tLN:1020778\n"
             << "@SQ\tSN:MU273392.1\tLN:189707\n"
             << "@SQ\tSN:ML143380.1\tLN:412368\n"
             << "@SQ\tSN:ML143378.1\tLN:461303\n"
             << "@SQ\tSN:KN196485.1\tLN:156562\n"
             << "@SQ\tSN:KQ458387.1\tLN:155930\n"
             << "@SQ\tSN:KQ458388.1\tLN:174749\n"
             << "@SQ\tSN:KN196486.1\tLN:153027\n"
             << "@SQ\tSN:KQ759761.1\tLN:145162\n"
             << "@SQ\tSN:ML143379.1\tLN:12295\n"
             << "@SQ\tSN:KQ759762.2\tLN:101040\n"
             << "@SQ\tSN:KV766199.1\tLN:188004\n"
             << "@SQ\tSN:ML143385.1\tLN:17435\n"
             << "@SQ\tSN:MU273397.1\tLN:330493\n"
             << "@SQ\tSN:ML143382.1\tLN:28824\n"
             << "@SQ\tSN:MU273393.1\tLN:68810\n"
             << "@SQ\tSN:ML143383.1\tLN:68192\n"
             << "@SQ\tSN:MU273394.1\tLN:140567\n"
             << "@SQ\tSN:ML143381.1\tLN:403128\n"
             << "@SQ\tSN:MU273396.1\tLN:294119\n"
             << "@SQ\tSN:ML143384.1\tLN:14678\n"
             << "@SQ\tSN:MU273395.1\tLN:619716\n"
             << "@SQ\tSN:KZ208923.1\tLN:48370\n"
             << "@SQ\tSN:MU273398.1\tLN:865743\n"
             << "@SQ\tSN:KZ208924.1\tLN:209722\n"
             << "@SQ\tSN:KN196487.1\tLN:101150\n"
             << "@SQ\tSN:KI270762.1\tLN:354444\n"
             << "@SQ\tSN:KI270766.1\tLN:256271\n"
             << "@SQ\tSN:KI270760.1\tLN:109528\n"
             << "@SQ\tSN:KI270765.1\tLN:185285\n"
             << "@SQ\tSN:GL383518.1\tLN:182439\n"
             << "@SQ\tSN:GL383519.1\tLN:110268\n"
             << "@SQ\tSN:GL383520.2\tLN:366580\n"
             << "@SQ\tSN:KI270764.1\tLN:50258\n"
             << "@SQ\tSN:KI270763.1\tLN:911658\n"
             << "@SQ\tSN:KI270759.1\tLN:425601\n"
             << "@SQ\tSN:KI270761.1\tLN:165834\n"
             << "@SQ\tSN:KI270770.1\tLN:136240\n"
             << "@SQ\tSN:KI270773.1\tLN:70887\n"
             << "@SQ\tSN:KI270774.1\tLN:223625\n"
             << "@SQ\tSN:KI270769.1\tLN:120616\n"
             << "@SQ\tSN:GL383521.1\tLN:143390\n"
             << "@SQ\tSN:KI270772.1\tLN:133041\n"
             << "@SQ\tSN:KI270775.1\tLN:138019\n"
             << "@SQ\tSN:KI270771.1\tLN:110395\n"
             << "@SQ\tSN:KI270768.1\tLN:110099\n"
             << "@SQ\tSN:GL582966.2\tLN:96131\n"
             << "@SQ\tSN:GL383522.1\tLN:123821\n"
             << "@SQ\tSN:KI270776.1\tLN:174166\n"
             << "@SQ\tSN:KI270767.1\tLN:161578\n"
             << "@SQ\tSN:JH636055.2\tLN:173151\n"
             << "@SQ\tSN:KI270783.1\tLN:109187\n"
             << "@SQ\tSN:KI270780.1\tLN:224108\n"
             << "@SQ\tSN:GL383526.1\tLN:180671\n"
             << "@SQ\tSN:KI270777.1\tLN:173649\n"
             << "@SQ\tSN:KI270778.1\tLN:248252\n"
             << "@SQ\tSN:KI270781.1\tLN:113034\n"
             << "@SQ\tSN:KI270779.1\tLN:205312\n"
             << "@SQ\tSN:KI270782.1\tLN:162429\n"
             << "@SQ\tSN:KI270784.1\tLN:184404\n"
             << "@SQ\tSN:KI270790.1\tLN:220246\n"
             << "@SQ\tSN:GL383528.1\tLN:376187\n"
             << "@SQ\tSN:KI270787.1\tLN:111943\n"
             << "@SQ\tSN:GL000257.2\tLN:586476\n"
             << "@SQ\tSN:KI270788.1\tLN:158965\n"
             << "@SQ\tSN:GL383527.1\tLN:164536\n"
             << "@SQ\tSN:KI270785.1\tLN:119912\n"
             << "@SQ\tSN:KI270789.1\tLN:205944\n"
             << "@SQ\tSN:KI270786.1\tLN:244096\n"
             << "@SQ\tSN:KI270793.1\tLN:126136\n"
             << "@SQ\tSN:KI270792.1\tLN:179043\n"
             << "@SQ\tSN:KI270791.1\tLN:195710\n"
             << "@SQ\tSN:GL383532.1\tLN:82728\n"
             << "@SQ\tSN:GL949742.1\tLN:226852\n"
             << "@SQ\tSN:KI270794.1\tLN:164558\n"
             << "@SQ\tSN:GL339449.2\tLN:1612928\n"
             << "@SQ\tSN:GL383530.1\tLN:101241\n"
             << "@SQ\tSN:KI270796.1\tLN:172708\n"
             << "@SQ\tSN:GL383531.1\tLN:173459\n"
             << "@SQ\tSN:KI270795.1\tLN:131892\n"
             << "@SQ\tSN:GL000250.2\tLN:4672374\n"
             << "@SQ\tSN:KI270800.1\tLN:175808\n"
             << "@SQ\tSN:KI270799.1\tLN:152148\n"
             << "@SQ\tSN:GL383533.1\tLN:124736\n"
             << "@SQ\tSN:KI270801.1\tLN:870480\n"
             << "@SQ\tSN:KI270802.1\tLN:75005\n"
             << "@SQ\tSN:KB021644.2\tLN:185823\n"
             << "@SQ\tSN:KI270797.1\tLN:197536\n"
             << "@SQ\tSN:KI270798.1\tLN:271782\n"
             << "@SQ\tSN:KI270804.1\tLN:157952\n"
             << "@SQ\tSN:KI270809.1\tLN:209586\n"
             << "@SQ\tSN:KI270806.1\tLN:158166\n"
             << "@SQ\tSN:GL383534.2\tLN:119183\n"
             << "@SQ\tSN:KI270803.1\tLN:1111570\n"
             << "@SQ\tSN:KI270808.1\tLN:271455\n"
             << "@SQ\tSN:KI270807.1\tLN:126434\n"
             << "@SQ\tSN:KI270805.1\tLN:209988\n"
             << "@SQ\tSN:KI270818.1\tLN:145606\n"
             << "@SQ\tSN:KI270812.1\tLN:282736\n"
             << "@SQ\tSN:KI270811.1\tLN:292436\n"
             << "@SQ\tSN:KI270821.1\tLN:985506\n"
             << "@SQ\tSN:KI270813.1\tLN:300230\n"
             << "@SQ\tSN:KI270822.1\tLN:624492\n"
             << "@SQ\tSN:KI270814.1\tLN:141812\n"
             << "@SQ\tSN:KI270810.1\tLN:374415\n"
             << "@SQ\tSN:KI270819.1\tLN:133535\n"
             << "@SQ\tSN:KI270820.1\tLN:36640\n"
             << "@SQ\tSN:KI270817.1\tLN:158983\n"
             << "@SQ\tSN:KI270816.1\tLN:305841\n"
             << "@SQ\tSN:KI270815.1\tLN:132244\n"
             << "@SQ\tSN:GL383539.1\tLN:162988\n"
             << "@SQ\tSN:GL383540.1\tLN:71551\n"
             << "@SQ\tSN:GL383541.1\tLN:171286\n"
             << "@SQ\tSN:GL383542.1\tLN:60032\n"
             << "@SQ\tSN:KI270823.1\tLN:439082\n"
             << "@SQ\tSN:GL383545.1\tLN:179254\n"
             << "@SQ\tSN:KI270824.1\tLN:181496\n"
             << "@SQ\tSN:GL383546.1\tLN:309802\n"
             << "@SQ\tSN:KI270825.1\tLN:188315\n"
             << "@SQ\tSN:KI270832.1\tLN:210133\n"
             << "@SQ\tSN:KI270830.1\tLN:177092\n"
             << "@SQ\tSN:KI270831.1\tLN:296895\n"
             << "@SQ\tSN:KI270829.1\tLN:204059\n"
             << "@SQ\tSN:GL383547.1\tLN:154407\n"
             << "@SQ\tSN:JH159136.1\tLN:200998\n"
             << "@SQ\tSN:JH159137.1\tLN:191409\n"
             << "@SQ\tSN:KI270827.1\tLN:67707\n"
             << "@SQ\tSN:KI270826.1\tLN:186169\n"
             << "@SQ\tSN:GL877875.1\tLN:167313\n"
             << "@SQ\tSN:GL877876.1\tLN:408271\n"
             << "@SQ\tSN:KI270837.1\tLN:40090\n"
             << "@SQ\tSN:GL383549.1\tLN:120804\n"
             << "@SQ\tSN:KI270835.1\tLN:238139\n"
             << "@SQ\tSN:GL383550.2\tLN:169178\n"
             << "@SQ\tSN:GL383552.1\tLN:138655\n"
             << "@SQ\tSN:GL383553.2\tLN:152874\n"
             << "@SQ\tSN:KI270834.1\tLN:119498\n"
             << "@SQ\tSN:GL383551.1\tLN:184319\n"
             << "@SQ\tSN:KI270833.1\tLN:76061\n"
             << "@SQ\tSN:KI270836.1\tLN:56134\n"
             << "@SQ\tSN:KI270840.1\tLN:191684\n"
             << "@SQ\tSN:KI270839.1\tLN:180306\n"
             << "@SQ\tSN:KI270843.1\tLN:103832\n"
             << "@SQ\tSN:KI270841.1\tLN:169134\n"
             << "@SQ\tSN:KI270838.1\tLN:306913\n"
             << "@SQ\tSN:KI270842.1\tLN:37287\n"
             << "@SQ\tSN:KI270844.1\tLN:322166\n"
             << "@SQ\tSN:KI270847.1\tLN:1511111\n"
             << "@SQ\tSN:KI270845.1\tLN:180703\n"
             << "@SQ\tSN:KI270846.1\tLN:1351393\n"
             << "@SQ\tSN:KI270852.1\tLN:478999\n"
             << "@SQ\tSN:KI270851.1\tLN:263054\n"
             << "@SQ\tSN:KI270848.1\tLN:327382\n"
             << "@SQ\tSN:GL383554.1\tLN:296527\n"
             << "@SQ\tSN:KI270849.1\tLN:244917\n"
             << "@SQ\tSN:GL383555.2\tLN:388773\n"
             << "@SQ\tSN:KI270850.1\tLN:430880\n"
             << "@SQ\tSN:KI270854.1\tLN:134193\n"
             << "@SQ\tSN:KI270856.1\tLN:63982\n"
             << "@SQ\tSN:KI270855.1\tLN:232857\n"
             << "@SQ\tSN:KI270853.1\tLN:2659700\n"
             << "@SQ\tSN:GL383556.1\tLN:192462\n"
             << "@SQ\tSN:GL383557.1\tLN:89672\n"
             << "@SQ\tSN:GL383563.3\tLN:375691\n"
             << "@SQ\tSN:KI270862.1\tLN:391357\n"
             << "@SQ\tSN:KI270861.1\tLN:196688\n"
             << "@SQ\tSN:KI270857.1\tLN:2877074\n"
             << "@SQ\tSN:JH159146.1\tLN:278131\n"
             << "@SQ\tSN:JH159147.1\tLN:70345\n"
             << "@SQ\tSN:GL383564.2\tLN:133151\n"
             << "@SQ\tSN:GL000258.2\tLN:1821992\n"
             << "@SQ\tSN:GL383565.1\tLN:223995\n"
             << "@SQ\tSN:KI270858.1\tLN:235827\n"
             << "@SQ\tSN:KI270859.1\tLN:108763\n"
             << "@SQ\tSN:GL383566.1\tLN:90219\n"
             << "@SQ\tSN:KI270860.1\tLN:178921\n"
             << "@SQ\tSN:KI270864.1\tLN:111737\n"
             << "@SQ\tSN:GL383567.1\tLN:289831\n"
             << "@SQ\tSN:GL383570.1\tLN:164789\n"
             << "@SQ\tSN:GL383571.1\tLN:198278\n"
             << "@SQ\tSN:GL383568.1\tLN:104552\n"
             << "@SQ\tSN:GL383569.1\tLN:167950\n"
             << "@SQ\tSN:GL383572.1\tLN:159547\n"
             << "@SQ\tSN:KI270863.1\tLN:167999\n"
             << "@SQ\tSN:KI270868.1\tLN:61734\n"
             << "@SQ\tSN:KI270865.1\tLN:52969\n"
             << "@SQ\tSN:GL383573.1\tLN:385657\n"
             << "@SQ\tSN:GL383575.2\tLN:170222\n"
             << "@SQ\tSN:GL383576.1\tLN:188024\n"
             << "@SQ\tSN:GL383574.1\tLN:155864\n"
             << "@SQ\tSN:KI270866.1\tLN:43156\n"
             << "@SQ\tSN:KI270867.1\tLN:233762\n"
             << "@SQ\tSN:GL949746.1\tLN:987716\n"
             << "@SQ\tSN:GL383577.2\tLN:128386\n"
             << "@SQ\tSN:KI270869.1\tLN:118774\n"
             << "@SQ\tSN:KI270871.1\tLN:58661\n"
             << "@SQ\tSN:KI270870.1\tLN:183433\n"
             << "@SQ\tSN:GL383578.2\tLN:63917\n"
             << "@SQ\tSN:KI270874.1\tLN:166743\n"
             << "@SQ\tSN:KI270873.1\tLN:143900\n"
             << "@SQ\tSN:GL383579.2\tLN:201197\n"
             << "@SQ\tSN:GL383580.2\tLN:74653\n"
             << "@SQ\tSN:GL383581.2\tLN:116689\n"
             << "@SQ\tSN:KI270872.1\tLN:82692\n"
             << "@SQ\tSN:KI270875.1\tLN:259914\n"
             << "@SQ\tSN:KI270878.1\tLN:186262\n"
             << "@SQ\tSN:KI270879.1\tLN:304135\n"
             << "@SQ\tSN:KI270876.1\tLN:263666\n"
             << "@SQ\tSN:KI270877.1\tLN:101331\n"
             << "@SQ\tSN:GL383583.2\tLN:96924\n"
             << "@SQ\tSN:GL383582.2\tLN:162811\n"
             << "@SQ\tSN:KI270880.1\tLN:284869\n"
             << "@SQ\tSN:KI270881.1\tLN:144206\n"
             << "@SQ\tSN:KI270892.1\tLN:162212\n"
             << "@SQ\tSN:KI270894.1\tLN:214158\n"
             << "@SQ\tSN:KI270893.1\tLN:161218\n"
             << "@SQ\tSN:KI270895.1\tLN:162896\n"
             << "@SQ\tSN:KI270896.1\tLN:378547\n"
             << "@SQ\tSN:KI270897.1\tLN:1144418\n"
             << "@SQ\tSN:KI270898.1\tLN:130957\n"
             << "@SQ\tSN:GL000251.2\tLN:4795265\n"
             << "@SQ\tSN:KI270899.1\tLN:190869\n"
             << "@SQ\tSN:KI270901.1\tLN:136959\n"
             << "@SQ\tSN:KI270900.1\tLN:318687\n"
             << "@SQ\tSN:KI270902.1\tLN:106711\n"
             << "@SQ\tSN:KI270903.1\tLN:214625\n"
             << "@SQ\tSN:KI270904.1\tLN:572349\n"
             << "@SQ\tSN:KI270906.1\tLN:196384\n"
             << "@SQ\tSN:KI270905.1\tLN:5161414\n"
             << "@SQ\tSN:KI270907.1\tLN:137721\n"
             << "@SQ\tSN:KI270910.1\tLN:157099\n"
             << "@SQ\tSN:KI270909.1\tLN:325800\n"
             << "@SQ\tSN:JH159148.1\tLN:88070\n"
             << "@SQ\tSN:KI270908.1\tLN:1423190\n"
             << "@SQ\tSN:KI270912.1\tLN:174061\n"
             << "@SQ\tSN:KI270911.1\tLN:157710\n"
             << "@SQ\tSN:GL949747.2\tLN:729520\n"
             << "@SQ\tSN:KB663609.1\tLN:74013\n"
             << "@SQ\tSN:KI270913.1\tLN:274009\n"
             << "@SQ\tSN:KI270924.1\tLN:166540\n"
             << "@SQ\tSN:KI270925.1\tLN:555799\n"
             << "@SQ\tSN:GL000252.2\tLN:4604811\n"
             << "@SQ\tSN:KI270926.1\tLN:229282\n"
             << "@SQ\tSN:KI270927.1\tLN:218612\n"
             << "@SQ\tSN:GL949748.2\tLN:1064304\n"
             << "@SQ\tSN:KI270928.1\tLN:176103\n"
             << "@SQ\tSN:KI270934.1\tLN:163458\n"
             << "@SQ\tSN:GL000253.2\tLN:4677643\n"
             << "@SQ\tSN:GL949749.2\tLN:1091841\n"
             << "@SQ\tSN:KI270935.1\tLN:197351\n"
             << "@SQ\tSN:GL000254.2\tLN:4827813\n"
             << "@SQ\tSN:GL949750.2\tLN:1066390\n"
             << "@SQ\tSN:KI270936.1\tLN:164170\n"
             << "@SQ\tSN:GL000255.2\tLN:4606388\n"
             << "@SQ\tSN:GL949751.2\tLN:1002683\n"
             << "@SQ\tSN:KI270937.1\tLN:165607\n"
             << "@SQ\tSN:GL000256.2\tLN:4929269\n"
             << "@SQ\tSN:GL949752.1\tLN:987100\n"
             << "@SQ\tSN:KI270758.1\tLN:76752\n"
             << "@SQ\tSN:GL949753.2\tLN:796479\n"
             << "@SQ\tSN:KI270938.1\tLN:1066800\n"
             << "@SQ\tSN:KI270882.1\tLN:248807\n"
             << "@SQ\tSN:KI270883.1\tLN:170399\n"
             << "@SQ\tSN:KI270884.1\tLN:157053\n"
             << "@SQ\tSN:KI270885.1\tLN:171027\n"
             << "@SQ\tSN:KI270886.1\tLN:204239\n"
             << "@SQ\tSN:KI270887.1\tLN:209512\n"
             << "@SQ\tSN:KI270888.1\tLN:155532\n"
             << "@SQ\tSN:KI270889.1\tLN:170698\n"
             << "@SQ\tSN:KI270890.1\tLN:184499\n"
             << "@SQ\tSN:KI270891.1\tLN:170680\n"
             << "@SQ\tSN:KI270914.1\tLN:205194\n"
             << "@SQ\tSN:KI270915.1\tLN:170665\n"
             << "@SQ\tSN:KI270916.1\tLN:184516\n"
             << "@SQ\tSN:KI270917.1\tLN:190932\n"
             << "@SQ\tSN:KI270918.1\tLN:123111\n"
             << "@SQ\tSN:KI270919.1\tLN:170701\n"
             << "@SQ\tSN:KI270920.1\tLN:198005\n"
             << "@SQ\tSN:KI270921.1\tLN:282224\n"
             << "@SQ\tSN:KI270922.1\tLN:187935\n"
             << "@SQ\tSN:KI270923.1\tLN:189352\n"
             << "@SQ\tSN:KI270929.1\tLN:186203\n"
             << "@SQ\tSN:KI270930.1\tLN:200773\n"
             << "@SQ\tSN:KI270931.1\tLN:170148\n"
             << "@SQ\tSN:KI270932.1\tLN:215732\n"
             << "@SQ\tSN:KI270933.1\tLN:170537\n"
             << "@SQ\tSN:GL000209.2\tLN:177381\n"
             << "@SQ\tSN:J01415.2\tLN:16569\n"
             << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188\tCL:bwa mem GCA_000001405.29_GRCh38.p14_genomic.fna SRR12922210_1.fastq\n"
             << "SRR12922210.1	0	CM000665.2	113658392	60	150M	*	0	0	CTACGACTACTGCCAGAGCTGGTTGACATGGGAGAGTCGGCCTGTCTACCGTTGATCAAAGAACCATTAATAACCTCGTGATCAGGGAGGCAGGAATGATGCTGTGGAGGATCTCTATCATCATTGTTCATCAGTAATAGTTCCTGTTTG	F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFFF:F,::FFFFFFFFF:FFFFFFFFFFFFFFF:FFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF	NM:i:1	MD:Z:77T72	AS:i:145	XS:i:0\n"
             << "SRR12922210.2	16	CM000669.2	129857657	60	150M	*	0	0	TCTAATAGATCTGTGAAAAAAAGAATCCCAATTGGGGAAAAGTATGATGAGGAACAGAATAGTACCAGGTCTTGAAGTGTCACCTTACAGACTGAAAGGGGAAAATAGTAACTATGCAATGGAGAAACCAGACAACACCTGAACCAGGTG	FFF,FFFFFFFFFFFFFFF:FFFFFFF,FFFF:FFFFFFFFFFFFFFFFFF,FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF	NM:i:0	MD:Z:150	AS:i:150	XS:i:0\n"
             << "SRR12922210.7	0	CM000664.2	201698711	60	150M	*	0	0	TGTAATCCACAGATGCCTAACCCAGCCTAACAACTAAATTACTCACAAGCCTTTCCTGCCACAAATCTTCTTAAGTATTCCTTGAAGTCTAATTTGTATCATTTTCAAACCTACATTCTAGTATTTTCTCTCAAATTGCCAACACATAGT	FF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:150	AS:i:150	XS:i:21\n"
             << "SRR12922210.8	16	CM000663.2	236409311	60	150M	*	0	0	TGTTTTGCTTTTTTTCTTTGTTATTTCCTTACGTTTTTATTACTCAGTATTTACATGTGCATCAGAAACCCCAAAGTAAGATTTTTCTAAGATAACTGCTTTATCTTGGTATTCTACGTGATAATAATAAACTTTGGGTTTCTCCAGGCT	FFFFFF:FFFFFFF:F:FFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F,F	NM:i:1	MD:Z:27T122	AS:i:145	XS:i:21\n"
             << "SRR12922210.9	0	CM000673.2	130462623	60	150M	*	0	0	TGCAGTACTGCACCAAGCTGTGGTGCACCGGGAAGGCCAAGGGACAGATGGTGTGCCAGACCCGCCACTTCCCCTGGGCCGATGGCACCAGCTGTGGCGAGGGCAAGCTCTGCCTCAAAGGGGCCTGCGTGGAGAGACACAACCTCAACA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF,FFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:0	MD:Z:150	AS:i:150	XS:i:0\n"
             << "SRR12922210.10	16	CM000672.2	100508473	60	150M	*	0	0	TAAGCACAGAAGAAAGAGAAAAGGCCAGGAAGCAAAAAGATGGGCCCAGTTACAGATTGACAGACAGCCTCTACAAGACAGAAGACCAGACCCCCAACTGCATGGCAGAGCCACATACAGGGAGAACAAAGCTGACCCCCCCCTCCATAA	FFFFFFFFFFFFFF:,FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFF,FFFFFFFFF,FFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF	NM:i:1	MD:Z:72C77	AS:i:145	XS:i:20\n"
             << "SRR12922210.1603	4	*	0	0	*	*	0	0	CACGGAGTTGTAGACTGAAATGCTGTCGATTGCGAAGGACATGCACATGAAGAAGACTAACATCTGCAGTGAAAACAAGAAGTAAGAGGCGTGAGGACCACTGGCTGCCGTAGAAGAGCAAGCAGGCTGACTGGGTACACAACCGTACAA	:FF,FF,FF,FF,F:,,:F:FFFFF,:FF,,,,,FFFFF,F,,,,,,,,FF,,F:,F,FFFFFFFFFF,F,FF:FF:F,:F,,FF:F:,:,:,F,FF,,,,,,FF,FF,FF,F::,F,F,,,,:,,,,,:F,,,,,F,F::F,,F:FFF,	AS:i:0	XS:i:0\n"
             << "SRR12922210.1867	4	*	0	0	*	*	0	0	TGCACCCCCAGAACCCCCCCCACCCCCGACCCCCGGGCCGGCCGCCCCCATGATTTTACCCACAGTTTAGAGATATTAAGGACACCACACAAAAAAACCAGTAAGAAAAGGTAAATTAATGGGACACGTACAGTTTGGGTTAACAAGGGG	FF,FFF,FFF,::F:FFF::F:FFFF,,,,,F:,,:F:F,,,F,:FFF,,,::,,FF,:F,,::,,,F,,F,,,:F,,,,F,,,,FF,,,FFF,F,,F,,:F,F,::F::,:,,:,F:,,,,,F,F,,,FFF,,::,FF,:,,,,,,FFF	AS:i:0	XS:i:0\n";
    }

    void loadSamData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        IOReader* pIoReader = new FileReader(filename);
        pIoReader->openIO(); 
        /*
         * 直接按原始内容载入整块（含 @SQ 头部），并记录换行符位置。
         * 头部在 reader 里独立成块后，actuator 测试需要一个自含头部的完整块。
         */
        const size_t blockSize = pInBlock->getBlockSize();
        const size_t readLen = pIoReader->readIO(pInBlock->getBuffer(), blockSize);
        pInBlock->setDataLen((int64_t)readLen);
        pInBlock->setBlockType(SAM);
        std::vector<size_t>& npos = pInBlock->getNpos();
        const char* buf = reinterpret_cast<const char*>(pInBlock->getBuffer());
        for (size_t i = 0; i < readLen; ++i) {
            if (buf[i] == '\n') {
                npos.push_back(i);
            }
        }
        delete(pIoReader);
    }

    void printBufferBinary(uint8_t* buffer, uint32_t bufferLen) {
        for (uint32_t j = 0; j < bufferLen; j++) {
            // Output src[i] in binary format
            for (int n = 7; n >= 0; n--) {
                fprintf(stderr, "%d", (buffer[j] >> n) & 1);
            }
            fprintf(stderr, "\t");

            if ((j+1) % 4 == 0) {
                fprintf(stderr, "\n");
            }
        }
        fprintf(stderr, "\n");
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
};

// Add actual test methods
// TEST_F(SamMappingTest, testMapingData) {
//     loadSamData(SamMappingData::testSamFile);
//     Reference refGene("../../data/GCA_000001405.29_GRCh38.p14_genomic.fna", 1);
//     refGene.makeIndex();
//     SamActuator actuator(pInBlock, pOutBlock, &refGene);
    
//     // Pre-analysis
//     int32_t result = actuator.preAnalysis();
//     EXPECT_EQ(result, 0);

//     uint32_t fieldSrcLen = 0;
//     Json::Value fieldMeta;

//     result = actuator.compressNumber<uint16_t>(1, fieldSrcLen, fieldMeta);
//     EXPECT_GT(result, 0);

//     result = actuator.compressChrName(2, fieldSrcLen, fieldMeta);
//     EXPECT_GT(result, 0);

//     result = actuator.compressNumber<int64_t>(3, fieldSrcLen, fieldMeta);
//     EXPECT_GT(result, 0);

//     uint32_t lineNum = pInBlock->npos.size();
//     uint8_t* buffer = pInBlock->getBuffer();
//     uint32_t fieldIdx = 9;
//     uint32_t offset = 0;
//     uint64_t nOffset = 0;

//     const uint32_t baseMaxLength = actuator.maxBaseLength + 4;
//     const uint32_t lsquash = (baseMaxLength >> 2) + !!(baseMaxLength & 0x3);
    
//     uint32_t baseMappedLength = (baseMaxLength << 1);
//     LOG_DEBUG("baseMappedLength = %d", baseMappedLength);
//     uint8_t* basePairBuffer = MemoryUtil::safeAlloc<uint8_t>(baseMaxLength);
//     uint8_t* baseSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(lsquash);
//     uint8_t* baseMappedBuffer = MemoryUtil::safeAlloc<uint8_t>(baseMappedLength);
//     actuator.baseNPosBuffer = MemoryUtil::safeAlloc<uint32_t>(actuator.baseNCount);


//     for (uint32_t lineIdx = actuator.headEndLine; lineIdx < lineNum; ++lineIdx) {
//         uint32_t lineStart = (lineIdx == 0) ? 0 : pInBlock->npos[lineIdx - 1] + 1;
//         uint32_t lineEnd = pInBlock->npos[lineIdx] - lineStart;
        
//         uint8_t* line = buffer + lineStart;
//         // Skip header lines (starting with @)
//         if (*line == '@') {
//             continue;
//         }
        
//         uint32_t contentIdx = lineIdx - actuator.headEndLine;
//         // Extract SEQ field (field 9)
//         if (fieldIdx >= actuator.contentPos[contentIdx].size() + 1) {
//             continue;
//         }
        
//         uint32_t prevTabPos = actuator.contentPos[contentIdx][fieldIdx - 1];
//         uint32_t currTabPos = (fieldIdx < actuator.contentPos[contentIdx].size()) ? actuator.contentPos[contentIdx][fieldIdx] : lineEnd;
//         uint8_t* seqStart = line + prevTabPos + 1;
//         uint32_t seqLength = currTabPos - prevTabPos - 1;
        
//         if (seqLength == 0)  {
//             continue;
//         }
        
//         // Get mapping information from SAM fields
//         uint16_t chrId = actuator.mappedChr[lineIdx];
//         uint64_t startPos = actuator.mappedPos[lineIdx];
        
//         // Extract FLAG field to determine strand
//         uint16_t flag = actuator.mappedFlag[lineIdx];
        
//         // Process sequence: remove N's and record positions
//         uint32_t outLen = 0;
//         uint32_t nCountInSeq = 0;
//         for (uint32_t n = 0; n < seqLength; n++) {
//             char ch = seqStart[n];
//             if (ch == 'N' || ch == 'n') {
//                 actuator.baseNPosBuffer[nOffset] = n;
//                 nOffset++;
//                 nCountInSeq++;
//             }
//         }
        
//         if (chrId != 0xFFFF && chrId != 0xFFFE) {
//             // Get chromosome start position from SamInfo
//             const ChromosomeInfo& chrInfo = SamInfo::getInstance().getChromosomeInfo(chrId);
//             uint64_t chrStartPos = chrInfo.position;
//             // Calculate actual reference position
//             int64_t refPos = chrStartPos + startPos - 1; // SAM is 1-based

//             uint8_t shiftBitLength = (refPos % 4);
//             int64_t refSquashPos = refPos / 4;
            
//             // Determine strand direction from FLAG bit 4
//             uint32_t squashBufferLength = 0;
//             bool isReverse = (flag & 0x10) != 0;
//             if (isReverse) {
//                 LOG_DEBUG("Revese: flag = %d", flag);
//                 // actgPair(basePairBuffer, seqStart, seqLength);
//                 squashBufferLength = actgSquash(seqStart, seqLength, baseSquashBuffer);
//             } else {
//                 squashBufferLength = actgSquash(seqStart, seqLength, baseSquashBuffer);
//             }

//             if (refGene.getSquashLength() < refSquashPos + squashBufferLength) {
//                 LOG_DEBUG("%ld, %ld", refGene.getSquashLength(), refSquashPos + squashBufferLength);
//                 continue;
//             }

//             uint8_t* refeMappedPos = MemoryUtil::safeAlloc<uint8_t>(squashBufferLength);
//             const uint8_t* beginRefPos = actuator.pRefeGene->getSquash() + refSquashPos;
//             if (shiftBitLength == 0) {
//                 memcpy(refeMappedPos, beginRefPos, squashBufferLength);
//             } else if (shiftBitLength == 1) {
//                 for (uint32_t i = 0; i < squashBufferLength; ++i) {
//                     refeMappedPos[i] = ((beginRefPos[i] << 2) & 0xFC) + ((beginRefPos[i + 1] >> 6) & 0x03);
//                 }
//             } else if (shiftBitLength == 2) {
//                 for (uint32_t i = 0; i < squashBufferLength; ++i) {
//                     refeMappedPos[i] = ((beginRefPos[i] << 4) & 0xF0) + ((beginRefPos[i + 1] >> 4) & 0x0F);
//                 }
//             } else if (shiftBitLength == 3) {
//                 for (uint32_t i = 0; i < squashBufferLength; ++i) {
//                     refeMappedPos[i] = ((beginRefPos[i] << 6) & 0xC0) + ((beginRefPos[i + 1] >> 2) & 0x3F);
//                 }
//             }

//             outLen = actgStretchMappingXor(baseSquashBuffer, refeMappedPos, squashBufferLength, baseMappedBuffer);
            
//             LOG_DEBUG("outlen = %d, squashBufferLength= %d, sequenLen = %d", outLen, squashBufferLength, seqLength);

//             LOG_DEBUG("baseSquashBuffer:");
//             printBufferBinary(baseSquashBuffer, squashBufferLength);

//             LOG_DEBUG("refeMappedPos:");
//             printBufferBinary(refeMappedPos, squashBufferLength);
            
//             LOG_DEBUG("baseMappedBuffer:");
//             printBufferBinary(baseMappedBuffer, seqLength);
            
//             MemoryUtil::safeFree(refeMappedPos);
//         } 
//         offset++;
//     }

//     MemoryUtil::safeFree(basePairBuffer);
//     MemoryUtil::safeFree(baseSquashBuffer);
//     MemoryUtil::safeFree(baseMappedBuffer);
//     MemoryUtil::safeFree(actuator.baseNPosBuffer);
// }

// TEST_F(SamMappingTest, testShirftMemory) {
//     uint8_t* src = MemoryUtil::safeAlloc<uint8_t>(16);
//     for (uint8_t i = 0; i < 16; ++i) {
//         src[i] = i;
//     }

//     for (int i = 0; i < 4; i++) {
//         for (int j = 0; j < 4; j++) {
//             // Output src[i] in binary format
//              for (int n = 7; n >= 0; n--) {
//                 printf("%d", (src[i*4 +j] >> n) & 1);
//             }
//             printf("\t");
//         }
//         printf("\n");
//     }

//     printf("\n");

//     uint8_t* dst = MemoryUtil::safeAlloc<uint8_t>(16);
//     for (int i = 0; i < 16; ++i) {
//         dst[i] = ((src[i] << 2) & 0xFC)  + ((src[(i + 1) % 16] >> 6) & 0x03 );
//     }

//     for (int i = 0; i < 4; i++) {
//         for (int j = 0; j < 4; j++) {
//             // Output src[i] in binary format
//              for (int n = 7; n >= 0; n--) {
//                 printf("%d", (dst[i*4 +j] >> n) & 1);
//             }
//             printf("\t");
//         }
//         printf("\n");
//     }

//     MemoryUtil::safeFree(src);
//     MemoryUtil::safeFree(dst);
// }