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

#include <Ptexture.h>
#include <PtexUtils.h>
#include <PtexHalf.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>
#include <time.h>
#include "KdTree.h"


class Vec3f
{
    float p[3];
    
 public:
    inline Vec3f()
    { p[0]=0.0f; p[1]=0.0f; p[2]=0.0f; }
    
    inline explicit Vec3f(float x)
    { p[0]=x; p[1]=x; p[2]=x; }
    
    inline Vec3f(float x, float y, float z)
    { p[0]=x; p[1]=y; p[2]=z; }
    
    inline Vec3f(const float *f)
    { p[0]=f[0]; p[1]=f[1]; p[2]=f[2]; }
    
    inline Vec3f(const double *d)
    { p[0]=d[0]; p[1]=d[1]; p[2]=d[2]; }
    
    inline float* getValue()
    { return p; }
    
    inline const float* getValue() const
    { return p; }
    
    inline float length() const
    { return std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]); }
    
    inline float& operator[](int n)
    { return p[n]; }
    
    inline const float& operator[](int n) const
    { return p[n]; }
    
    inline Vec3f operator+(const Vec3f &v) const
    { return Vec3f(p[0]+v.p[0], p[1]+v.p[1], p[2]+v.p[2]); }
    
    inline Vec3f operator-(const Vec3f &v) const
    { return Vec3f(p[0]-v.p[0], p[1]-v.p[1], p[2]-v.p[2]); }
    
    inline Vec3f operator-() const
    { return Vec3f(-p[0], -p[1], -p[2]); }
    
    inline Vec3f operator*(float f) const
    { return Vec3f(p[0]*f, p[1]*f, p[2]*f); } const
    
    inline Vec3f operator/(float f)
    { return Vec3f(p[0]/f, p[1]/f, p[2]/f); }
};


class PointCloud
{
 private:
    unsigned char* _values;
    int _totalpoints;
    int _datasize;
    KdTree<3> _tree;
        
 public:
    PointCloud(int datasize) : _values(0), _totalpoints(0), _datasize(datasize) {}
    ~PointCloud() { if (_values) delete [] _values; }
    
    bool setSize(int numpts)
    {
        _totalpoints = numpts;
        if (_values) delete [] _values;
        _values = new unsigned char[_totalpoints*_datasize];
        _tree.allocatePoints(_totalpoints);
        return true;
    }
    
    int size() { return _totalpoints; }
    
    void sortTree() { _tree.sort(); }
    
    void store(int i, const float position[3], const void* data)
    {
        if (i>=_totalpoints) return;
        _tree.setPoint(i, (float*)&position[0]);
        memcpy(&_values[i*_datasize], data, _datasize);
    }
    
    bool find(const float position[3], void* data, float distance = 1e9f)
    {
        int index = _tree.findNearest(&position[0], distance);
        if (index<0) return false;
        memcpy(data, &_values[_tree.id(index)*_datasize], _datasize);
        return true;
    }
};


class Mesh
{
    std::vector<Vec3f> verts;
    std::vector< std::vector<int> > faceVertIndices;

 public:
    
    int numFaces() const { return faceVertIndices.size(); }
    
    inline bool evaluate(int faceid, float u, float v, Vec3f& p, int resU = 0, int resV = 0)
    {
        if (faceid>=int(faceVertIndices.size())) return false;
        
        std::vector<int>& vind = faceVertIndices[faceid];
        
        if (vind.size()==4)
        {
            p = verts[vind[0]]*(1-u)*(1-v) +
                verts[vind[1]]*u*(1-v) +
                verts[vind[2]]*u*v +
                verts[vind[3]]*(1-u)*v;
            return true;
        }
        
        if (vind.size()==3)
        {
            if (resU>0 && resV>0)
            {
                if (u+v<=1)
                {
                    u = u - 1.0f/(resU*6);
                    v = v - 1.0f/(resV*6);
                }
                else
                {
                    u = 1 + 1.0f/(resU*6) - u;
                    v = 1 + 1.0f/(resV*6) - v;
                }
            }
            p = verts[vind[0]] +
                (verts[vind[1]]-verts[vind[0]])*u + 
                (verts[vind[2]]-verts[vind[0]])*v;
            return true;
        }
        
        return false;
    }
    
    float faceTexelSize(int faceid, int resU, int resV)
    {
        std::vector<int>& vind = faceVertIndices[faceid];
        if (vind.size()==4)
        {
            return std::min(
                std::min((verts[vind[1]] - verts[vind[0]]).length(),
                         (verts[vind[2]] - verts[vind[3]]).length()) / resU,
                std::min((verts[vind[3]] - verts[vind[0]]).length(),
                         (verts[vind[2]] - verts[vind[1]]).length()) / resV);
                
        }
        else if (vind.size()==3)
        {
            return std::min( std::min((verts[vind[1]] - verts[vind[0]]).length(),
                         (verts[vind[2]] - verts[vind[1]]).length()),
                         (verts[vind[0]] - verts[vind[2]]).length()) / resU;
        }
        
        return 0;
    }
    
    Ptex::Res calculateFaceResolution(int faceid, float texelSize)
    {
        std::vector<int>& vind = faceVertIndices[faceid];
        if (vind.size()==4)
        {
            uint32_t resU = PtexUtils::floor_log2(std::max((verts[vind[1]] - verts[vind[0]]).length(),
                            (verts[vind[2]] - verts[vind[3]]).length()) / texelSize);
            uint32_t resV = PtexUtils::floor_log2(std::max((verts[vind[3]] - verts[vind[0]]).length(),
                            (verts[vind[2]] - verts[vind[1]]).length()) / texelSize);
            
            return Ptex::Res(resU,resV);
                
        }
        else if (vind.size()==3)
        {
            float len = std::max(std::max((verts[vind[1]] - verts[vind[0]]).length(),
                                 (verts[vind[2]] - verts[vind[1]]).length()),
                                 (verts[vind[0]] - verts[vind[2]]).length());
            
            uint32_t resU = PtexUtils::floor_log2(len / texelSize);
            return Ptex::Res(resU,resU);
        }
        
        return Ptex::Res(1,1);
    }
    
    bool getMeshData(PtexTexture* r)
    {
        if (!r) return false;

        PtexPtr<PtexMetaData> meta(r->getMetaData());
        if (!meta) return false;

        if (meta->numKeys()<3)
        {
            std::cerr << "Geometry meta data not found in Ptex file\n";
            return false;
        }

        const float* vp;
        const int *fvi, *fvc;
        int vertcount, indexcount, facevertcount;

        meta->getValue("PtexFaceVertCounts",fvc,facevertcount);
        meta->getValue("PtexVertPositions",vp,vertcount);
        meta->getValue("PtexFaceVertIndices",fvi,indexcount);

        if (facevertcount==0)
        {
            std::cerr << "Missing meta data : PtexFaceVertCounts\n";
            return false;
        }

        if (vertcount==0)
        {
            std::cerr << "Missing meta data : PtexVertPositions\n";
            return false;
        }

        if (indexcount==0)
        {
            std::cerr << "Missing meta data : PtexFaceVertIndices\n";
            return false;
        }

        verts.resize(vertcount/3);
        for (uint i=0; i<verts.size(); i++) verts[i] = Vec3f(vp[i*3],vp[i*3+1],vp[i*3+2]);

        int vi = 0;
        for (int fi=0; fi<facevertcount; fi++)
        {
            int nverts = fvc[fi];
            if (r->meshType()==Ptex::mt_quad && nverts!=4)
            {
                // Split n-sided polygon into quads

                Vec3f cp = Vec3f(0,0,0); // center of polygon
                int vib = verts.size();

                // Find edge midpoints and polygon center
                for (int i=0; i<nverts; i++)
                {
                    Vec3f& v = verts[fvi[vi+i]];
                    cp = cp + v;
                    verts.push_back((v+verts[fvi[vi+(i+1)%nverts]])/2);
                }
                verts.push_back(cp/nverts);

                // Add new quads using edge midpoints and polygon center
                for (int i=0; i<nverts; i++)
                {
                    std::vector<int> v;
                    v.push_back(fvi[vi+i]);
                    v.push_back(vib+i);
                    v.push_back(verts.size()-1);
                    v.push_back(vib+(i+nverts-1)%nverts);
                    faceVertIndices.push_back(v);
                }

                vi += nverts;
            }
            else if (r->meshType()==Ptex::mt_triangle && nverts!=3)
            {
                std::cerr << "Illegal vertex count for triangular Ptex\n";
                return false;
            }
            else
            {
                std::vector<int> v;
                for (int i=0; i<nverts; i++) v.push_back(fvi[vi++]);
                faceVertIndices.push_back(v);
            }
        }

        if (int(faceVertIndices.size())!=r->numFaces())
        {
            std::cerr << "Face count mismatch\n";
            return false;
        }

        return true;
    }
};


class PixelBuffer
{
    int _numChannels;
    Ptex::DataType _dataType;
    void* _data;
    
 public:
    PixelBuffer(Ptex::DataType dataType, int numChannels, int numPixels) : _numChannels(numChannels), _dataType(dataType), _data(0)
    {
        int numValues = numPixels*_numChannels;
        
        switch (_dataType)
        {
            case Ptex::dt_uint8:  _data = new uint8_t[numValues];  break;
            case Ptex::dt_float:  _data = new float[numValues];    break;
            case Ptex::dt_half:   _data = new PtexHalf[numValues]; break;
            case Ptex::dt_uint16: _data = new uint16_t[numValues]; break;
        }
    }
    
    ~PixelBuffer()
    {
        switch (_dataType)
        {
            case Ptex::dt_uint8:  delete [] ((uint8_t*)_data);  break;
            case Ptex::dt_float:  delete [] ((float*)_data);    break;
            case Ptex::dt_half:   delete [] ((PtexHalf*)_data); break;
            case Ptex::dt_uint16: delete [] ((uint16_t*)_data); break;
        }
    }

    inline void writePixel(float* val, int pixelIndex)
    {
        int d = _numChannels*pixelIndex;
        
        switch (_dataType)
        {
        case Ptex::dt_uint8:
            for (int i=0; i<_numChannels; i++) ((uint8_t*)_data)[d+i] = uint8_t(val[i]*255);
            break;
        case Ptex::dt_float:
            for (int i=0; i<_numChannels; i++) ((float*)_data)[d+i] = val[i];
            break;
        case Ptex::dt_half:
            for (int i=0; i<_numChannels; i++) ((PtexHalf*)_data)[d+i] = val[i];
            break;
        case Ptex::dt_uint16:
            for (int i=0; i<_numChannels; i++) ((uint16_t*)_data)[d+i] = uint16_t(val[i]*65535);
        }
    }
    
    void* getData() { return _data; }
    
    void* getPixel(int pixelIndex)
    { return (unsigned char*)_data+pixelIndex*_numChannels*Ptex::DataSize(_dataType); }
};


void imageMirrorRows(unsigned char* data, int w, int h, int datasize)
{
    unsigned char tmp[datasize];
    int halfw = w/2;
    
    for (int pix=0; pix<h*halfw; pix++)
    {
        int y = pix/halfw;
        int x = pix-y*halfw;
        unsigned char* row = data+y*w*datasize;
        memcpy(tmp, row+x*datasize, datasize);
        memcpy(row+x*datasize, row+(w-1-x)*datasize, datasize);
        memcpy(row+(w-1-x)*datasize, tmp, datasize);
    }
}

void imageRotate(unsigned char* data, int w, int h, int datasize, int rotations)
{
    if (rotations!=1 && rotations!=3) return;
    
    unsigned char* newdata = new unsigned char[w*h*datasize];
    
    if (rotations==1)
    {
        for (int y=0; y<h; y++)
            for (int x=0; x<w; x++)
                memcpy(newdata+(y+(w-1-x)*h)*datasize, data+(x+y*w)*datasize, datasize);
    }
    else if (rotations==3)
    {
        for (int y=0; y<h; y++)
            for (int x=0; x<w; x++)
                memcpy(newdata+((h-1-y)+x*h)*datasize, data+(x+y*w)*datasize, datasize);
    }
    memcpy(data, newdata, datasize*w*h);
    delete [] newdata;
}

bool transferPtex(std::string input, std::string output, float searchDist = -1, float matchDist = -1)
{
    time_t timer1;
    time(&timer1);
    
    std::string ptexError;
    
    PtexPtr<PtexTexture> fromPtex(PtexTexture::open(input.c_str(), ptexError, true));
    if (!fromPtex)
    {
        std::cerr << ptexError << "\n";
        return false;
    }
    
    PtexPtr<PtexTexture> toPtex(PtexTexture::open(output.c_str(), ptexError, true));
    if (!toPtex)
    {
        std::cerr << ptexError << "\n";
        std::cerr << "Output file must exist and contain geometry meta data for transfer\n";
        return false;
    }
    
    if (fromPtex->meshType()!=toPtex->meshType())
    {
        std::cerr << "Mesh types are inconsistent (quad / triangle)\n";
        return false;
    }
    
    Mesh fromMesh;
    Mesh toMesh;
    if (!fromMesh.getMeshData(fromPtex) || !toMesh.getMeshData(toPtex))
    {
        return false;
    }
    
    Ptex::MeshType meshType = fromPtex->meshType();
    
    Ptex::DataType dataType = fromPtex->dataType();
    int targetAlpha = fromPtex->alphaChannel();
    
    int numchan = fromPtex->numChannels();
    int datasize = numchan*Ptex::DataSize(fromPtex->dataType());
    PtexMetaData* meta = toPtex->getMetaData();
    
    // Find total number of texels
    int totalTexels = 0;
    std::vector<int> baseIndex(fromPtex->numFaces());
    for (int faceid=0; faceid<fromPtex->numFaces(); faceid++)
    {
        baseIndex[faceid] = totalTexels;
        totalTexels += fromPtex->getFaceInfo(faceid).res.size();
    }
    float avgTexelSize = 0;
    
    // Create point cloud
    PointCloud fullPts(datasize);
    fullPts.setSize(totalTexels);
    
    // Sparse point cloud with four samples for each face, for face matching
    PointCloud facePts(sizeof(int)*2);
    facePts.setSize(fromPtex->numFaces()*4);
    
    // Populate point cloud using input Ptex data
    for (int faceid=0; faceid<fromPtex->numFaces(); faceid++)
    {
        std::cout << "  Reading faces: " << faceid+1 << " / " << fromPtex->numFaces() << "\r" << std::flush;
        
        const Ptex::FaceInfo& f = fromPtex->getFaceInfo(faceid);
        
        unsigned char* data = new unsigned char[f.res.size()*datasize];
        fromPtex->getData(faceid, data, 0);
        
        // Add points for each texel in input face
#pragma omp parallel for
        for (int pix=0; pix<f.res.size(); pix++)
        {
            int uPixel = pix%f.res.u();
            int vPixel = pix/f.res.u();
            float u = (float(uPixel)+0.5f)/f.res.u();
            float v = (float(vPixel)+0.5f)/f.res.v();
            Vec3f p;
            fromMesh.evaluate(faceid,u,v,p,f.res.u(),f.res.v());
            fullPts.store(baseIndex[faceid] + pix, p.getValue(), data+pix*datasize);
        }
        delete [] data;
        
        float s = fromMesh.faceTexelSize(faceid, f.res.u(), f.res.v());
        avgTexelSize += s;
    }
    
    if (toPtex->meshType()==Ptex::mt_quad)
    {
#pragma omp parallel for
        for (int faceid=0; faceid<fromPtex->numFaces(); faceid++)
        {
            // Add four samples for each input face to sparse point cloud
            Vec3f p;
            int faceData[] = { faceid,0 };
            for (int i=0; i<4; i++)
            {
                float u = (i==0 || i==3) ? 0.25 : 0.75;
                float v = (i<2) ? 0.25 : 0.75;
                faceData[1] = i;
                fromMesh.evaluate(faceid, u, v, p, 2, 2);
                facePts.store(faceid*4+i, p.getValue(), &faceData);
            }
        }
    }
    
    avgTexelSize /= fromPtex->numFaces();
    
    std::cout << "\n  Sorting texels...\n";
    
    fullPts.sortTree();
    if (toPtex->meshType()==Ptex::mt_quad) facePts.sortTree();
    
    if (searchDist<0)
    {
        searchDist = avgTexelSize*40;
        std::cout << "  Calculating a search distance: " << searchDist << "\n";
    }
    if (matchDist<0) matchDist = avgTexelSize;
    
    // Write output data
    
    PtexPtr<PtexWriter> w(PtexWriter::open(output.c_str(), meshType, dataType,
                          numchan, targetAlpha, toMesh.numFaces(), ptexError));

    float errval[numchan];
    for (int i=0; i<numchan; i++) errval[i] = 1;
    
    int copied = 0;
        
    for (int faceid=0; faceid<toMesh.numFaces(); faceid++)
    {
        std::cout << "  Writing faces: " << faceid+1 << " / " << toMesh.numFaces() << "\r" << std::flush;
        
        Ptex::FaceInfo f = toPtex->getFaceInfo(faceid);
        
        // Check for direct face copy
        int fromFace = -1;
        int faceData[4][2];
        
        if (matchDist>0 && toPtex->meshType()==Ptex::mt_quad)
        {
            Vec3f pos;
            for (int i=0; i<4; i++)
            {
                float u = (i==0 || i==3) ? 0.25 : 0.75;
                float v = (i<2) ? 0.25 : 0.75;
                if (!toMesh.evaluate(faceid,u,v,pos,2,2)) { fromFace = -1; break; }
                if (!facePts.find(pos.getValue(), &faceData[i], matchDist)) { fromFace = -1; break; }
                if (i==0) fromFace = faceData[i][0];
                else if (faceData[i][0]!=fromFace) { fromFace = -1; break; }
            }
        }
        
        if (fromFace>=0)
        {
            bool match = true;
            bool flipU = false;
            bool flipV = false;
            int rotations = 0;
            
            const int& v0 = faceData[0][1];
            const int& v1 = faceData[1][1];
            const int& v2 = faceData[2][1];
            const int& v3 = faceData[3][1];
            
            if (v0==0 && v1==1 && v2==2 && v3==3)
            {
                // Direct match, no flip/rotate needed
            }
            else if (v0==3 && v1==2 && v2==1 && v3==0)
            {
                flipV = true;
            }
            else if (v0==1 && v1==0 && v2==3 && v3==2)
            {
                flipU = true;
            }
            else if (v0==2 && v1==3 && v2==0 && v3==1)
            {
                flipU = flipV = true;
            }
            else if (v0==1 && v1==2 && v2==3 && v3==0)
            {
                rotations = 1;
            }
            else if (v0==3 && v1==0 && v2==1 && v3==2)
            {
                rotations = 3;
            }
            else if (v0==2 && v1==1 && v2==0 && v3==3)
            {
                flipV = true;
                rotations = 1;
            }
            else if (v0==0 && v1==3 && v2==2 && v3==1)
            {
                flipV = true;
                rotations = 3;
            }
            else
            {
                match = false;
            }
            
            if (match)
            {
                f.res = fromPtex->getFaceInfo(fromFace).res;
                unsigned char* data = new unsigned char[f.res.size()*datasize];
                
                if (flipV)
                    fromPtex->getData(fromFace, data+f.res.u()*datasize*(f.res.v()-1), -f.res.u()*datasize);
                else
                    fromPtex->getData(fromFace, data, 0);
                
                if (flipU) imageMirrorRows(data, f.res.u(), f.res.v(), datasize);
                if (rotations==1||rotations==3)
                {
                    imageRotate(data, f.res.u(), f.res.v(), datasize, rotations);
                    f.res.swapuv();
                }
                
                w->writeFace(faceid, f, data, 0);
                delete [] data;
                
                copied++;

                continue;
            }
        }
        
        f.res = toMesh.calculateFaceResolution(faceid,avgTexelSize);
        
        PixelBuffer pixels(dataType,numchan,f.res.size());
        
#pragma omp parallel for
        for (int pix=0; pix<f.res.size(); pix++)
        {
            float u = (float(pix%f.res.u())+0.5f)/f.res.u();
            float v = (float(pix/f.res.u())+0.5f)/f.res.v();
            Vec3f p;
            if (searchDist>0 && toMesh.evaluate(faceid,u,v,p,f.res.u(),f.res.v()))
            {
                if (fullPts.find(p.getValue(), pixels.getPixel(pix), searchDist)) continue;
            }
            pixels.writePixel(&errval[0], pix);
        }
        
        w->writeFace(faceid, f, pixels.getData(), 0);
    }
    
    std::cout << "\n";
    if (copied>0) std::cout << "  Copied " << copied << " faces using a tolerance of " << matchDist << "\n";
    
    std::cout << "  Adding meta data...\n";
    w->writeMeta(meta);
    
    std::cout << "  Closing file...\n";
    if (!w->close(ptexError))
    {
        std::cerr << ptexError << std::endl;
        return false;
    }
    
    time_t timer2;
    time(&timer2);
    double seconds = difftime(timer2,timer1);
    if (seconds<2)
        std::cout << "  Finished";
    else
        std::cout << "  Finished in " << seconds << " seconds";
        
    int threads = omp_get_max_threads();
    std::cout << " using " << threads << " thread" << ((threads>1)?"s":"") << "\n";
    
    return true;
}
   
void showUsage()
{
    std::cout << "Usage:\n";
    std::cout << "    ptxtransfer [options] <old.ptx> <new.ptx>\n\n";
    std::cout << "Options:\n";
    std::cout << "    -d distance    maximum texel search distance\n";
    std::cout << "    -m distance    tolerance for direct face transfer\n";
    std::cout << "    -h             show help message\n";
    std::cout << "    -t num         set number of threads\n";
}

int main(int argc, char **argv)
{
    if (argc<3)
    {
        showUsage();
        return 1;
    }
    
    std::string input = "";
    std::string output = "";
    int numthreads = 1;
    float distance = -1;
    float match = -1;
    
    for (int i = 1; i < argc;)
    {
        std::string arg = argv[i++];
        
        if (arg[0] == '-')
        {
            if (arg=="-h" || arg=="-help")
            {
                showUsage();
                return 0;
            }
            else if (arg=="-t")
            {
                if (i>=argc) { showUsage(); return 1; }
                numthreads = atoi(argv[i++]);
            }
            else if (arg=="-d")
            {
                if (i>=argc) { showUsage(); return 1; }
                distance = atof(argv[i++]);
            }
            else if (arg=="-m")
            {
                if (i>=argc) { showUsage(); return 1; }
                match = atof(argv[i++]);
            }
            else
            {
                std::cerr << "Error: Invalid option: " << arg << "\n";
                return 1;
            }
        }
        else if (input=="") input = arg;
        else if (output=="") output = arg;
        else 
        {
            showUsage();
            return 0;
        }
    }
    
    if (input=="" || output=="")
    {
        showUsage();
        return 0;
    }
    
    omp_set_num_threads(numthreads);
    transferPtex(input,output,distance,match);
}
