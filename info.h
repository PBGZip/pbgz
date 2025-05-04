#ifndef _INFO_H_
#define _INFO_H_

#include <iostream>
#if __x86_64__  || __i386__ || _M_X64 || _M_IX86
#include "hardware.h"
#endif

typedef enum
{
    ZIP,
    UNZIP,
    MODE_UNKNOW
} mode;

class info
{
public:
    info()
    {
        reset();
    }

    virtual ~info()
    {
    }

    void reset()
    {
        refgene.clear();
        m = MODE_UNKNOW;
        force = false;
        decompres2gz = false;
        nopackref = false;
        remove = false;
        threads = 0;
        clevel = 0;
        cmd.clear();
#if __x86_64__  || __i386__ || _M_X64 || _M_IX86
        support_simd = hardware().support_simd();
#else
        support_simd = false;
#endif
    }

    info& operator=(const info& data)
    {
        m = data.m;
        infile = data.infile;
        outfile = data.outfile;
        refgene = data.refgene;
        remove = data.remove;
        nopackref = data.nopackref;
        force = data.force;
        decompres2gz = data.decompres2gz;
        threads = data.threads;
        clevel = data.clevel;
        support_simd = data.support_simd;
        cmd = data.cmd;
        return *this;
    }

    void set_remove(const bool &r) {this->remove = r;}

    bool get_remove() const {return this->remove;}
    
    void set_clevel(const int &l) {this->clevel = l;}

    int get_clevel() const {return this->clevel;}

    void set_threads(const int &t) {this->threads = t;}

    int get_threads() const {return this->threads;}

    void set_nopackref(const bool &n) {this->nopackref = n;}

    bool get_nopackref() const {return this->nopackref;}

    void set_refgene(const std::string &r) {this->refgene = r;}

    const std::string &get_refgene() const {return this->refgene;}
    
    void set_force(const bool &f) {this->force = f;}

    bool get_force() const {return this->force;}

    void set_decompres2gz(const bool &to_gz) {this->decompres2gz = to_gz;}

    bool get_decompress2gz() const { return this->decompres2gz; }

    void set_mode(const mode &m) {this->m = m;}

    const mode & get_mode() const {return this->m;}

    void set_infile(const std::string &input) {this->infile = input;}

    const std::string & get_infile() const {return this->infile;}

    void set_outfile(const std::string &output) {this->outfile = output;}

    const std::string & get_outfile() const {return this->outfile;}

    bool get_support_simd() const { return this->support_simd; }

    void set_cmd(const std::string cmd) {this->cmd = cmd;}

    const std::string &get_cmd() const {return this->cmd;}

    void dump()  /* debug */
    {
        return;
        fprintf(stderr, "mode: %s\n", (m == ZIP) ? "ZIP" : ((m == UNZIP) ? "UNZIP" : "MODE_UNKNOW"));
        fprintf(stderr, "cmd: %s\n", cmd.c_str());
        fprintf(stderr, "infile: %s\n", infile.c_str());
        fprintf(stderr, "outfile: %s\n", outfile.c_str());
        fprintf(stderr, "refgene: %s\n", refgene.c_str());
        fprintf(stderr, "remove: %d\n", remove);
        fprintf(stderr, "nopackref: %d\n", nopackref);
        fprintf(stderr, "force: %d\n", force);
        fprintf(stderr, "remove: %d\n", remove);
        fprintf(stderr, "decompres2gz: %d\n", decompres2gz);
        fprintf(stderr, "threads: %d\n", threads);
        fprintf(stderr, "clevel: %d\n\n", clevel);
        fprintf(stderr, "support_simd: %d\n\n", support_simd);
    }

private:
    mode m;
    std::string infile, outfile;
    std::string refgene;
    std::string cmd;

    bool remove; // reomve input file
    bool nopackref;
    bool force;
    bool decompres2gz;
    int threads;
    int clevel; // compress level

    bool support_simd; // 支持硬件加速指令
};

#endif