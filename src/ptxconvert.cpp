/*
* Copyright Disney Enterprises, Inc.  All rights reserved.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License
* and the following modification to it: Section 6 Trademarks.
* deleted and replaced with:
*
* 6. Trademarks. This License does not grant permission to use the
* trade names, trademarks, service marks, or product names of the
* Licensor and its affiliates, except as required for reproducing
* the content of the NOTICE file.
*
* You may obtain a copy of the License at
* http://www.apache.org/licenses/LICENSE-2.0
*/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include "Ptexture.h"
#include "PtexUtils.h"
#include "PtexHalf.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#if __GNUC__ >= 9
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#include <OpenImageIO/imageio.h>


struct Options {
    Ptex::BorderMode uMode, vMode;
    Ptex::DataType dt;
    std::string envfaces[6];
    bool src_is_envcube, dst_is_envcube;
    std::string src, dst;
    std::map<std::string,std::string> meta;

    Options() : 
        uMode(Ptex::m_clamp), vMode(Ptex::m_clamp),
        dt(Ptex::DataType(0)),
        src_is_envcube(0), dst_is_envcube(0) {}
} opt;


OIIO_NAMESPACE_USING

template<typename T>
inline void divalpha(T* data, int npixels, int nchannels, int alphachan, float scale)
{
    int alphaoffset; // offset to alpha chan from data ptr
    int nchandiv;    // number of channels to alpha-divide
    if (alphachan == 0) {
        // first channel is alpha chan: div the rest of the channels
        data++;
        alphaoffset = -1;
        nchandiv = nchannels - 1; 
    }
    else {
        // div all channels up to alpha chan
        alphaoffset = alphachan;
        nchandiv = alphachan;
    }

    for (T* end = data + npixels*nchannels; data != end; data += nchannels) {
        T alpha = data[alphaoffset];
        if (!alpha) continue; // don't divide by zero!
        float aval = scale / alpha;
        if (scale>1) aval = std::min(aval,1.0f); // clamp Int8 or Int16 data
        for (int i = 0; i < nchandiv; i++)      data[i] = T(data[i] * aval);
    }
}


void divAlpha(void* data, int npixels, Ptex::DataType dt, int nchannels, int alphachan)
{
    float scale = Ptex::OneValue(dt);
    switch(dt) {
    case Ptex::dt_uint8:  ::divalpha((uint8_t*) data, npixels, nchannels, alphachan, scale); break;
    case Ptex::dt_uint16: ::divalpha((uint16_t*) data, npixels, nchannels, alphachan, scale); break;
    case Ptex::dt_half:   ::divalpha((PtexHalf*) data, npixels, nchannels, alphachan, scale); break;
    case Ptex::dt_float:  ::divalpha((float*) data, npixels, nchannels, alphachan, scale); break;
    }
}

static const float uint8Scale = 1.0/UINT8_MAX;
static const float uint16Scale = 1.0/UINT16_MAX;

struct Img {
    // Simple class to pass raw image data
    // Data is assumed to be bottom-up and interleaved with no row padding
    Ptex::DataType dt;
    int nchan, achan;
    int w, h;
    void* data;
    std::string path;
    
    Img() : dt(Ptex::DataType(0)), nchan(0), achan(0),
              w(0), h(0), data(0) {}
    ~Img() { if (data) free(data); }
    
    bool load(const std::string& filename)
    {
        std::unique_ptr<ImageInput> in(ImageInput::create(filename));
        if (!in) return false;
        ImageSpec spec;
        if (!in->open(filename, spec))
        {
            std::cerr << "Unable to open file " << filename << std::endl;
            in->close();
            return 0;
        }
        
        if (spec.format==TypeDesc::UINT8) dt = Ptex::dt_uint8;
        else if (spec.format==TypeDesc::UINT16) dt = Ptex::dt_uint16;
        else if (spec.format==TypeDesc::FLOAT) dt = Ptex::dt_float;
        else if (spec.format==TypeDesc::HALF) dt = Ptex::dt_half;
        else
        {
            std::cerr << "Unsupported image data type" << std::endl;
            in->close();
            return 0;
        }
        
        nchan = spec.nchannels;
        achan = spec.alpha_channel;
        w = spec.width;
        h = spec.height;
        path = filename;
        
        int size = w * h * nchan * Ptex::DataSize(dt);
        data = malloc(size);
        
        int ystride = spec.nchannels*spec.format.size()*spec.width;
        
        if (!in->read_image(spec.format, ((char*)(data))+ystride*(h-1), AutoStride, -ystride))
        {
            std::cerr << "Unable to read image data" << std::endl;
            in->close();
            return 0;
        }

        if (achan != -1) divAlpha(data, w*h, dt, nchan, achan);
        
        in->close();
        
        return true;
    }
    
    bool save(const std::string& filename)
    {
        std::unique_ptr<ImageOutput> out(ImageOutput::create(filename));
        if (!out) return false;
        
        TypeDesc type = TypeDesc::UINT8;
        switch (dt)
        {
        case Ptex::dt_uint8:  type = TypeDesc::UINT8; break;
        case Ptex::dt_uint16: type = TypeDesc::UINT16; break;
        case Ptex::dt_float:  type = TypeDesc::FLOAT; break;
        case Ptex::dt_half:  type = TypeDesc::HALF; break;
        }
        ImageSpec spec(w, h, nchan, type);
        
        if (!out->open(filename, spec))
        {
            std::cerr << "Unable to save file " << filename << std::endl;
            out->close();
            return 0;
        }
        
        int ystride = spec.nchannels*spec.format.size()*spec.width;
        
        if (!out->write_image(spec.format, ((char*)(data))+ystride*(h-1), AutoStride, -ystride))
        {
            std::cerr << "Unable to write image data" << std::endl;
            out->close();
            return 0;
        }
        
        return true;
    }
    
    void toFloat(int u, int v, float* p)
    {
        int d = (u+v*w)*nchan;
        switch (dt)
        {
        case Ptex::dt_uint8:
            for (int i=0; i<nchan; i++) p[i] = ((uint8_t*)data)[d+i]*uint8Scale;
            break;
        case Ptex::dt_uint16:
            for (int i=0; i<nchan; i++) p[i] = ((uint16_t*)data)[d+i]*uint16Scale;
            break;
        case Ptex::dt_float:
            for (int i=0; i<nchan; i++) p[i] = ((float*)data)[d+i];
            break;
        case Ptex::dt_half:
            for (int i=0; i<nchan; i++) p[i] = ((PtexHalf*)data)[d+i];
            break;
        }
        
        if (achan != -1)
        {
            float a = p[achan];
            for (int i=0; i<nchan; i++) if (i!=achan) p[i] *= a;
        }
    }
    
};

bool PtexToImg(PtexTexture* ptx, Img& img, int faceid, bool flip)
{
    if (faceid < 0 || faceid >= ptx->numFaces()) {
        std::cerr << "Invalid faceid: " << faceid << std::endl;
        return 0;
    }

    img.dt = ptx->dataType();
    img.nchan = ptx->numChannels();
    img.achan = ptx->alphaChannel();
    Ptex::FaceInfo fi = ptx->getFaceInfo(faceid);
    img.w = fi.res.u();
    img.h = fi.res.v();
    int rowlen = img.w * img.nchan * Ptex::DataSize(img.dt);
    int size = rowlen * img.h;
    img.data = malloc(size);
    char* data = (char*) img.data;
    int stride = rowlen;
    if (flip) {
        data += rowlen * (img.h-1);
        stride = -rowlen;
    }
    ptx->getData(faceid, data, stride);
    return 1;
}


void usage()
{
    std::cerr << "Usage:\n"
              << "    ptxconvert [options] srcfile dstfile\n"
              << "\n"
              << "Options:\n"
              << "    -mode (clamp|black|periodic)\n"
              << "    -umode (clamp|black|periodic)\n"
              << "    -vmode (clamp|black|periodic)\n"
              << "    -envcube pxfile nxfile pyfile nyfile pzfile nzfile\n"
              << "    -meta key value\n"
              << "\n"
              << "Note:\n"
              << "    -envcube can be used in place of srcfile or dstfile\n"
              << "\n"
              << "Supported file types:\n"
              << "    .tif, .jpg, .png, etc. - standard image files\n"
              << "    .tx - prman (tif) texture\n"
              << "    .ptx - ptex texture (single face only)\n"
              << "    .penv - ptex environment cube map\n"
              << "\n"
              << "Not yet supported:\n"
              << "    .env - prman (tif) environment cube map\n"
              << "\n";
}


void usagehint()
{
    std::cerr << "(See ptxconvert -h for usage info)\n";
}


void error(const std::string& err)
{
    std::cerr << "Error: " << err << "\n";
    usagehint();
}


void badoption(const std::string& opt)
{
    std::cerr << "Error: Invalid option: " << opt << "\n";
    usagehint();
}


void missingoptval(const std::string& opt)
{
    std::cerr << "Error: Missing value for option: " << opt << "\n";
    usagehint();
}


bool checkPowerOfTwo(const Img& img)
{
    if (!PtexUtils::isPowerOfTwo(img.w) || !PtexUtils::isPowerOfTwo(img.h)) {
        std::cerr << "Error: Image resolution not power-of-two: "
                  << img.w << "x" << img.h << " " << img.path << "\n";
        return 0;
    }
    return 1;
}


bool endsWith(const std::string& str, const std::string& end)
{
    int strlen = str.length();
    int endlen = end.length();
    return strlen >= endlen && str.substr(strlen-endlen) == end;
}


bool readImage(const std::string& path, Img& img)
{
    if (endsWith(path, ".ptx")) {
        img.path = path;
        // read ptx face-zero image
        Ptex::String error;
        PtexTexture* ptx = PtexTexture::open(path.c_str(), error);
        if (!ptx) {
            std::cerr << error << std::endl;
            return 0;
        }
        if (ptx->numFaces() != 1) {
            std::cerr << "Error: Ptex file has multiple faces (not yet supported)" << std::endl;
            ptx->release();
            return 0;
        }
        bool ok = PtexToImg(ptx, img, /*faceid*/ 0, /*flip*/ true);
        ptx->release();
        return ok;
    }
    else {
        return img.load(path);
    }
}


bool writeImage(const std::string& path, Img& img)
{
    if (endsWith(path, ".ptx")) {
        // make sure it's a power of two
        if (!checkPowerOfTwo(img)) return 0;

        // write ptx face-zero image
        Ptex::String error;
        PtexWriter* ptx = PtexWriter::open(path.c_str(), Ptex::mt_quad, img.dt,
                                           img.nchan, img.achan, 1, error);
        if (!ptx) {
            std::cerr << error << std::endl;
            return 0;
        }
        ptx->setBorderModes(opt.uMode, opt.vMode);
        // write image top-down (flipped)
        int rowlen = img.w * img.nchan * Ptex::DataSize(img.dt);
        char* toprow = (char*) img.data + (img.h-1) * rowlen;
        Ptex::Res res(PtexUtils::floor_log2(img.w), PtexUtils::floor_log2(img.h));
        int adjfaces[4], adjedges[4];
        if (opt.uMode == Ptex::m_periodic) { adjfaces[0] = 0; adjfaces[2] = 0; adjedges[0] = 2; adjedges[2] = 0; }
        else { adjfaces[0] = -1; adjfaces[2] = -1; adjedges[0] = 0; adjedges[2] = 0; }
        if (opt.vMode == Ptex::m_periodic) { adjfaces[1] = 0; adjfaces[3] = 0; adjedges[1] = 3; adjedges[3] = 1; }
        else { adjfaces[1] = -1; adjfaces[3] = -1; adjedges[1] = 0; adjedges[3] = 0; }
        Ptex::FaceInfo finfo(res, adjfaces, adjedges);
        ptx->writeFace(0, finfo, toprow, -rowlen);

        // write meta data
        for (std::map<std::string,std::string>::iterator iter = opt.meta.begin(), end = opt.meta.end();
             iter != end; iter++)
        {
            ptx->writeMeta(iter->first.c_str(), iter->second.c_str());
        }

        bool ok = ptx->close(error);
        if (!ok) {
            std::cerr << error << std::endl;
        }
        ptx->release();
        return ok;
    }
    else {
        return img.save(path);
    }
}


bool readEnvcube(const std::string& path, Img images[6])
{
    Ptex::String error;
    PtexTexture* ptx = PtexTexture::open(path.c_str(), error);
    if (!ptx) {
        std::cerr << error << std::endl;
        return 0;
    }
    if (ptx->numFaces() != 6) {
        std::cerr << "Error: Ptex envcube file doesn't have six faces: "<< path << "\n";
        ptx->release();
        return 0;
    }
    bool ok;
    for (int i = 0; i < 6; i++) {
        ok = PtexToImg(ptx, images[i], /*faceid*/ i, /*flip*/ false);
        if (!ok) break;
    }
    ptx->release();
    return ok;
}


bool writeEnvcube(const std::string& path, Img images[6])
{
    if (endsWith(path, ".env")) {
        std::cerr << "Writing prman envcube files (.env) not yet supported: " << path << "\n";
        return 0;
    }
    if (!endsWith(path, ".penv")) {
        std::cerr << "Unsupported envcube file format: " << path << "\n";
        return 0;
    }

    // check that all faces have the same format, and power-of-two size
    Ptex::DataType dt = images[0].dt;
    int nchan = images[0].nchan;
    int achan = images[0].achan;
    for (int i = 1; i < 6; i++) {
        if (images[i].dt != dt || images[i].nchan != nchan || images[i].achan != achan) {
            std::cerr << "Cube face image format mismatch\n";
            return 0;
        }
        if (!checkPowerOfTwo(images[i])) return 0;
    }

    // build face info
    int adjfaces[6][4] = {
        { 3, 5, 2, 4 }, // px
        { 3, 4, 2, 5 }, // nx
        { 4, 0, 5, 1 }, // py
        { 5, 0, 4, 1 }, // ny
        { 3, 0, 2, 1 }, // pz
        { 3, 1, 2, 0 }  // nz
    };
    int adjedges[6][4] = {
        { 1, 3, 1, 1 }, // px
        { 3, 3, 3, 1 }, // nx
        { 2, 2, 2, 2 }, // py
        { 0, 0, 0, 0 }, // ny
        { 2, 3, 0, 1 }, // pz
        { 0, 3, 2, 1 }  // nz
    };
    Ptex::FaceInfo fi[6];
    for (int i = 0; i < 6; i++) {
        int w = images[i].w, h = images[i].h;
        Ptex::Res res(PtexUtils::floor_log2(w), PtexUtils::floor_log2(h));
        fi[i] = Ptex::FaceInfo(res, adjfaces[i], adjedges[i]);
    }

    // write file
    Ptex::String error;
    PtexWriter* w = PtexWriter::open(opt.dst.c_str(), Ptex::mt_quad, dt, nchan, achan,
                                     6, error);
    if (!w) {
        std::cerr << error << std::endl;
        return 0;
    }
    for (int i = 0; i < 6; i++) {
        w->writeFace(i, fi[i], images[i].data);
    }

    bool ok = w->close(error);
    if (!ok) std::cerr << error << "\n";
    w->release();

    return ok;
}


bool parseArgs(int argc, char** argv)
{
    // parse args
    for (int i = 1; i < argc;) {
        std::string arg = argv[i++], val;
        if (arg[0] == '-') {
            if (arg.length() == 1) { usage(); return 0; }
            if (arg == "-h" || arg == "-help") {
                usage();
                return 0;
            }
            else if (arg == "-e" || arg == "-envcube") {
                if (argc - i < 6) { missingoptval(arg); return 0; }
                for (int n = 0; n < 6; n++) opt.envfaces[n] = argv[i+n];
                i += 6;
                if (opt.src.empty()) {
                    opt.src_is_envcube = 1;
                    opt.src = "envcube";
                }
                else if (opt.dst.empty()) {
                    opt.dst_is_envcube = 1;
                    opt.dst = "envcube";
                }
                else {
                    error("-envcube must be in place of srcfile or dstfile");
                    return 0;
                }
                if (opt.src_is_envcube && opt.dst_is_envcube) {
                    error("-envcube may not be used for both srcfile and dstfile");
                    return 0;
                }
            }
            else if (arg == "-meta") {
                if (argc - i < 2) { missingoptval(arg); return 0; }
                opt.meta[argv[i]] = argv[i+1];
                i += 2;
            }
            else if (arg == "-mode" || arg == "-umode" || arg == "-vmode") {
                std::string modestr = argv[i++];
                Ptex::BorderMode mode;
                if (modestr == "clamp") mode = Ptex::m_clamp;
                else if (modestr == "black") mode = Ptex::m_black;
                else if (modestr == "periodic") mode = Ptex::m_periodic;
                else {
                    std::string err = "Unrecognized border mode value, '";
                    err += modestr;
                    err += "', for option ";
                    err += arg;
                    error(err);
                    return 0;
                }
                if (arg == "-umode" || arg == "-mode") opt.uMode = mode;
                if (arg == "-vmode" || arg == "-mode") opt.vMode = mode;
            }
            else {
                badoption(arg);
                return 0;
            }
        }
        else {
            if (opt.src.empty()) opt.src = arg;
            else if (opt.dst.empty()) opt.dst = arg;
            else {
                error("Too many arguments");
                return 0;
            }
        }
    }
    if (opt.src.empty()) { error("Missing srcfile argument"); return 0; }
    if (opt.dst.empty()) { error("Missing dstfile argument"); return 0; }
    return 1;
}


bool buildEnvcube()
{
    if (!endsWith(opt.dst, ".penv")) {
        std::cerr << "Unsupported envmap format: " << opt.dst << std::endl;
        return 0;
    }

    Img images[6];
    for (int i = 0; i < 6; i++) {
        if (!readImage(opt.envfaces[i], images[i])) return 0;
    }

    if (!writeEnvcube(opt.dst, images)) return 0;
    return 1;
}


bool extractEnvcube()
{
    // read envmap
    Img images[6];
    if (!readEnvcube(opt.src, images)) return 0;

    // write individual image files
    for (int i = 0; i < 6; i++) {
        if (!writeImage(opt.envfaces[i], images[i]))
            return 0;
    }
    return 1;
}


bool copyEnvcube()
{
    // read envmap
    Img images[6];
    if (!readEnvcube(opt.src, images)) return 0;

    // write envmap
    if (!writeEnvcube(opt.dst, images)) return 0;
    return 1;
}


bool convertFile(const std::string& src, const std::string& dst)
{
    if (endsWith(opt.src, ".penv") && endsWith(opt.dst, ".penv")) {
        return copyEnvcube();
    }

    Img img;
    if (!readImage(src, img)) return 0;
    if (!writeImage(dst, img)) return 0;
    return 1;
}


int main(int argc, char** argv)
{
    if (!parseArgs(argc, argv)) return 1;

    bool ok = 1;
    if (opt.src_is_envcube) {
        if (opt.dst_is_envcube) {
            std::cerr << "Can't use -envcube for both src and dst" << std::endl;
            ok = 0;
        }
        else
            ok = buildEnvcube();
    }
    else if (opt.dst_is_envcube) {
        ok = extractEnvcube();
    }
    else {
        ok = convertFile(opt.src, opt.dst);
    }
    return ok ? 0 : 1;
}
