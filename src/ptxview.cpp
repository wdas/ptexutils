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


#define GL_GLEXT_PROTOTYPES

#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <GL/freeglut.h>
#else
# include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif

#ifndef GL_HALF_FLOAT_ARB
#define GL_HALF_FLOAT_ARB 0x140B
#endif

#include "Ptexture.h"

#define MENU_BGCOLOR 1
#define MENU_GRIDLINES 2
#define MENU_FOCUSVIEW 3
#define MENU_RESETVIEW 4
#define MENU_FULLSCREEN 5
#define MENU_QUIT 6
#define INVALIDNAME 9999999
#define DEGRAD (180.0/M_PI)

class Vec3f
{
    float p[3];
    
public:
    Vec3f() { p[0]=0.0; p[1]=0.0; p[2]=0.0; }
    Vec3f(float x, float y, float z) { p[0]=x; p[1]=y; p[2]=z; }
    Vec3f(const float *f) { p[0]=f[0]; p[1]=f[1]; p[2]=f[2]; }
    inline float* getValue() { return p; }
    inline double length() { return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]); }
    float dot(const Vec3f &v) const
    { return p[0]*v.p[0] + p[1]*v.p[1] + p[2]*v.p[2]; }
    float& operator[](int n)
    { return p[n]; }
    const float& operator[](int n) const
    { return p[n]; }
    Vec3f operator+(const Vec3f &v)
    { return Vec3f(p[0]+v.p[0], p[1]+v.p[1], p[2]+v.p[2]); }
    Vec3f operator-(const Vec3f &v)
    { return Vec3f(p[0]-v.p[0], p[1]-v.p[1], p[2]-v.p[2]); }
    Vec3f operator-()
    { return Vec3f(-p[0], -p[1], -p[2]); }
    Vec3f operator*(double f)
    { return Vec3f(p[0]*f, p[1]*f, p[2]*f); }
    Vec3f operator/(double f)
    { return Vec3f(p[0]/f, p[1]/f, p[2]/f); }
    Vec3f& operator+=(const Vec3f& v)
    { p[0]+=v.p[0]; p[1]+=v.p[1]; p[2]+=v.p[2]; return *this; }
    Vec3f& operator*=(double f)
    { p[0]*=f; p[1]*=f; p[2]*=f; return *this; }
    Vec3f& operator/=(double f)
    { p[0]/=f; p[1]/=f; p[2]/=f; return *this; }
};

struct BBox3f
{
    Vec3f bmin, bmax;
    
    BBox3f(const Vec3f& a, const Vec3f& b)
    {
        bmin = a;
        bmax = b;
    }
    
    BBox3f() {}

    inline Vec3f& getBound(int i) { return (i==0) ? bmin : bmax; }
    
    inline void grow(Vec3f& pt)
    {
        for (int i=0; i<3; i++)
        {
            bmin[i] = std::min(pt[i], bmin[i]);
            bmax[i] = std::max(pt[i], bmax[i]);
        }
    }
    
    inline void grow(BBox3f& b)
    {
        for (int i=0; i<3; i++)
        {
            bmin[i] = std::min(b.bmin[i], bmin[i]);
            bmax[i] = std::max(b.bmax[i], bmax[i]);
        }
    }
    
    inline Vec3f getCenter() { return (bmin+bmax)/2.0; }
    inline float diagonal() { return (bmax-bmin).length(); }
};


class camera
{
    double dist, azim, elev;
    Vec3f lookat;
    
public:
    camera();
    void applyProjectionTransform(float aspect, const BBox3f& box);
    void applyModelViewTransform();
    void rotateXY(float x, float y);
    void pan(float x, float y, bool oriented = true);
    void zoom(int z);
    float getDistance() { return dist; }
    Vec3f getLookAt() { return lookat; }
    void setLookAt(Vec3f v) { lookat = v; }
    void setDistance(float d) { dist = d; }
    void reset() { azim = elev = 0.0; lookat = Vec3f(0,0,0); dist = 10.0; };
};

camera::camera() : dist(10.0), azim(0), elev(0)
{
    lookat = Vec3f(0,0,0);
}

void camera::applyProjectionTransform(float aspect, const BBox3f& box)
{
    Vec3f lookdir;
    lookdir[0] = sin(azim/DEGRAD);
    lookdir[2] = -cos(azim/DEGRAD);
    lookdir[1] = sin(elev/DEGRAD);
    lookdir[0] *= cos(elev/DEGRAD);
    lookdir[2] *= cos(elev/DEGRAD);
    Vec3f lookpos = lookat-(lookdir*dist);
    float nearClip = 10e8;
    float farClip = -10e8;

    for (int i=0; i<8; i++)
    {
        Vec3f corner = Vec3f(((i/4)%2)?box.bmin[0]:box.bmax[0],
                             ((i/2)%2)?box.bmin[1]:box.bmax[1],
                             (i%2)?box.bmin[2]:box.bmax[2]);
        float depth = lookdir.dot(corner - lookpos);
        nearClip = std::min(nearClip, depth);
        farClip = std::max(farClip, depth);
    }
    nearClip -= nearClip/1000; farClip += farClip/1000;
    nearClip = std::max(nearClip, std::max(float(farClip/1000),1e-2f));

    double fov = 2.0*atan((18.0/std::min(aspect,1.0f))/35.0)*DEGRAD;
    gluPerspective(fov,aspect,nearClip,farClip);
}

void camera::applyModelViewTransform()
{
    glTranslated(0.0, 0.0, -dist);
    glRotated(-elev, 1.0, 0.0, 0.0);
    glRotated(azim, 0.0, 1.0, 0.0);
    glTranslated(-lookat[0],-lookat[1],-lookat[2]);
}

void camera::rotateXY(float x, float y)
{
    azim += x;
    elev += y;
}

void camera::pan(float x, float y, bool oriented)
{
    if (oriented)
    {
        lookat[0] -= cos(azim/DEGRAD)*x;
        lookat[2] -= sin(azim/DEGRAD)*x;
        lookat[1] -= cos(elev/DEGRAD)*y;
        lookat[0] -= sin(elev/DEGRAD)*y*sin(-azim/DEGRAD);
        lookat[2] -= sin(elev/DEGRAD)*y*cos(azim/DEGRAD);
    }
    else lookat += Vec3f(-x,-y,0);
}

void camera::zoom(int z)
{
    if (z<0) dist*=1.0+(((float)(-z))/50.0);
    else dist/=1.0+(((float)(z))/50.0);
    dist = std::max(dist,1e-4);
    dist = std::min(dist,1e5);
}

class mesh;

struct texturePage
{
    GLuint id;
    std::vector<int> faceIDs;
    
    texturePage();
    void deleteData();
    void makeGLTexture(mesh* m, PtexTexture* r, GLenum format, GLenum typegl);
};

texturePage::texturePage()
{
    id = 0;
}

void texturePage::deleteData()
{
    if (glIsTexture(id)) glDeleteTextures(1, &id);
}

struct polygon
{
    int indices[4];
    unsigned short ures, vres, ub, vb;
    void drawFace(std::vector<Vec3f>& verts, float id);
    void drawFlat(float id, float asp);
};

void polygon::drawFace(std::vector<Vec3f>& verts, float id)
{
    glTexCoord3f(0,0,id);
    glVertex3fv(verts[indices[0]].getValue());
    glTexCoord3f(1,0,id);
    glVertex3fv(verts[indices[1]].getValue());
    if (indices[3]>=0) glTexCoord3f(1,1,id);
    else glTexCoord3f(0,1,id);
    glVertex3fv(verts[indices[2]].getValue());
    if (indices[3]<0) return;
    glTexCoord3f(0,1,id);
    glVertex3fv(verts[indices[3]].getValue());
}

void polygon::drawFlat(float id, float asp)
{
    float xb = 0, yb = 0;
    float xd = 1, yd = 1;
    
    if (indices[3]>=0)
    {
        if (asp>1)
        {
            yb = 0.5 - 0.5/asp;
            yd = 1/asp;
        }
        else if (asp<1)
        {
            xb = 0.5 - 0.5*asp;
            xd = asp;
        }
    }

    glTexCoord3f(0,0,id);
    glVertex3f(xb,yb,0);
    glTexCoord3f(1,0,id);
    glVertex3f(xb+xd,yb,0);
    if (indices[3]<0)
    {
        glTexCoord3f(0,1,id);
        glVertex3f(0.5,0.866025404,0);
        return;
    }
    glTexCoord3f(1,1,id);
    glVertex3f(xb+xd,yb+yd,0);
    glTexCoord3f(0,1,id);
    glVertex3f(xb,yb+yd,0);
}

class mesh
{
    std::vector<Vec3f> verts;
    std::vector<polygon> polys;
    std::vector<texturePage> pages;
    bool triangles;
   
public:
    mesh() { triangles = false; }
    void makePages(PtexTexture* r);
    void draw(int faceID = -1);
    void deleteData();
    bool loadFromMetaData(PtexTexture* r);
    void makeFlatGeom(PtexTexture* r);
    bool empty() { return verts.size()==0 || polys.size()==0; }
    polygon* getFace(int fi)
        { return (fi>=0 && fi<(int)polys.size()) ? &polys[fi] : 0; }
    BBox3f getBounds();
    void addFace(polygon& q) { polys.push_back(q); }
    void addVertex(const float* f) { verts.push_back(Vec3f(f[0],f[1],f[2])); }
    inline Vec3f facePosition(int fi)
        { return (verts[polys[fi].indices[0]]+verts[polys[fi].indices[2]])/2; }
};

void mesh::makeFlatGeom(PtexTexture* r)
{
    int numfaces = r->numFaces();
    triangles = r->meshType()==Ptex::mt_triangle;
    
    float x=0, y=0;
    
    verts.clear();
    polys.clear();
    
    if (triangles)
    {
        int rowlen = (int)sqrt(numfaces);
        for (int i=0; i<numfaces; i++)
        {
            verts.push_back(Vec3f(x,y,0));
            verts.push_back(Vec3f(x+1,y,0));
            verts.push_back(Vec3f(x+0.5,y+0.866025404,0));

            polygon q;
            q.indices[0] = verts.size()-3;
            q.indices[1] = verts.size()-2;
            q.indices[2] = verts.size()-1;
            q.indices[3] = -1;
            addFace(q);

            x+=1/0.866025404;
            if (x>rowlen/0.866025404)
            {
                x = 0;
                y += 1;
            }
        }
    }
    else
    {
        int rs = (numfaces<2) ? 2 : (int)ceil(sqrt(numfaces))+1;
        for (int i=0; i<numfaces; i++)
        {
            x = i%(rs-1);
            y = i/(rs-1);

            polygon q;
            q.indices[0] = y*rs+x;
            q.indices[1] = y*rs+x+1;
            q.indices[2] = (y+1)*rs+x+1;
            q.indices[3] = (y+1)*rs+x;
            addFace(q);
            
            if ((int)verts.size()<=q.indices[2]) verts.resize(q.indices[2]+1);
            
            verts[q.indices[0]] = Vec3f(x,y,0);
            verts[q.indices[1]] = Vec3f(x+1,y,0);
            verts[q.indices[2]] = Vec3f(x+1,y+1,0);
            verts[q.indices[3]] = Vec3f(x,y+1,0);
        }
    }
}

bool mesh::loadFromMetaData(PtexTexture* r)
{
    PtexMetaData* meta = r->getMetaData();
    if (meta->numKeys()<3) return false;
    
    triangles = r->meshType()==Ptex::mt_triangle;

    const float* vf;
    const int *vif, *fvc;
    int vertcount, indexcount, facevertcount;

    meta->getValue("PtexFaceVertCounts",fvc,facevertcount);
    if (facevertcount==0) return false;
    meta->getValue("PtexVertPositions",vf,vertcount);
    if (vertcount==0) return false;
    meta->getValue("PtexFaceVertIndices",vif,indexcount);
    if (indexcount==0) return false;
    
    for (int i=0; i<vertcount; i+=3) addVertex(&vf[i]);
        
    int* fviptr = (int*) vif;
    int* nvptr = (int*) fvc;
    for (int i=0; i<facevertcount; i++)
    {
        int nverts = *nvptr;
        if (triangles)
        {
            if (nverts!=3) return false;
            polygon q;
            q.indices[0] = fviptr[0];
            q.indices[1] = fviptr[1];
            q.indices[2] = fviptr[2];
            q.indices[3] = -1;
            addFace(q);
        }
        else if (nverts==4)
        {
            polygon q;
            q.indices[0] = fviptr[0];
            q.indices[1] = fviptr[1];
            q.indices[2] = fviptr[2];
            q.indices[3] = fviptr[3];
            addFace(q);
        }
        else
        {
            int vib = verts.size();

            // Add new vertices at edge midpoints and polygon center
            Vec3f cp = Vec3f(0,0,0);
            for (int v=0; v<nverts; v++)
            {
                verts.push_back((verts[fviptr[v]] + verts[fviptr[(v+1)%nverts]])/2);
                cp += verts[fviptr[v]];
            }
            verts.push_back(cp/nverts);

            // Add the new sub-faces
            for (int v=0; v<nverts; v++)
            {
                polygon q;
                q.indices[0] = fviptr[v];
                q.indices[1] = vib+v;
                q.indices[2] = verts.size()-1;
                q.indices[3] = vib+(v+nverts-1)%nverts;
                addFace(q);
            }
        }
        fviptr+=nverts;
        nvptr++;
    }
    
    if ((int)polys.size()!=r->numFaces())
    {
        verts.clear();
        polys.clear();
        return false;
    }
    
    return true;
}

void mesh::makePages(PtexTexture* r)
{
    std::map< std::pair<int,int>, std::vector<int> > sizeMap;
    
    for (unsigned int fi=0; fi<polys.size(); fi++)
    {
        std::pair<int,int> res(polys[fi].ures, polys[fi].vres);
        sizeMap[res].push_back(fi);
    }
    
    GLint maxFaceCount;
    glGetIntegerv(0x88FF,&maxFaceCount);
    
    pages.clear();
    
    std::map< std::pair<int,int>, std::vector<int> >::iterator it = sizeMap.begin();
    while (it!=sizeMap.end())
    {
        std::vector<int>& fl = it->second;
        
        texturePage pg;
        for (unsigned int i=0; i<fl.size(); i++)
        {
            if ((int)pg.faceIDs.size()>=maxFaceCount)
            {
                pages.push_back(pg);
                pg.faceIDs.clear();
            }
            pg.faceIDs.push_back(fl[i]);
        }
        pages.push_back(pg);
        
        it++;
    }
    
    GLenum format, typegl;
    
    switch (r->dataType())
    {
        case Ptex::dt_uint16:
            typegl = GL_UNSIGNED_SHORT;
            break;
        case Ptex::dt_float:
            typegl = GL_FLOAT;
            break;
        case Ptex::dt_half:
            typegl = GL_HALF_FLOAT_ARB;
            break;
        default:
            typegl = GL_UNSIGNED_BYTE;
            break;
    }
    switch (r->numChannels())
    {
        case 1: format = GL_LUMINANCE; break;
        case 2: format = GL_LUMINANCE_ALPHA; break;
        case 3: format = GL_RGB; break;
        case 4: format = GL_RGBA; break;
        default: format = GL_LUMINANCE; break;
    }
    
    for (unsigned int p = 0; p<pages.size(); p++)
        pages[p].makeGLTexture(this, r, format, typegl);
}

void mesh::draw(int faceID)
{
    GLint rm;
    glGetIntegerv(GL_RENDER_MODE, &rm);
    bool selectMode = (rm==GL_SELECT);
    
    for (unsigned int p=0; p<pages.size(); p++)
    {
        texturePage& pg = pages[p];
        
        if (selectMode)
        {
            for (unsigned int i=0; i<pg.faceIDs.size(); i++)
            {
                glLoadName(pg.faceIDs[i]);
                glBegin(triangles?GL_TRIANGLES:GL_QUADS);
                polys[pg.faceIDs[i]].drawFace(verts,i);
                glEnd();
            }
            continue;
        }
        
        glBindTexture(GL_TEXTURE_2D_ARRAY,pg.id);
        glBegin(triangles?GL_TRIANGLES:GL_QUADS);
        for (unsigned int i=0; i<pg.faceIDs.size(); i++)
        {
            int f = pg.faceIDs[i];
            if (faceID>=0 && f!=faceID) continue;
            
            if (faceID>=0||polys.size()==1)
                polys[f].drawFlat(i,(float)polys[f].ures/(float)polys[f].vres);
            else
                polys[f].drawFace(verts,i);
        }
        glEnd();
    }
    
    if (selectMode) glLoadName(INVALIDNAME);
}

void mesh::deleteData()
{
    for (unsigned int p=0; p<pages.size(); p++) pages[p].deleteData();
}

BBox3f mesh::getBounds()
{
    BBox3f b(Vec3f(1e8,1e8,1e8),Vec3f(-1e8,-1e8,-1e8));
    
    for (unsigned int v=0; v<verts.size(); v++) b.grow(verts[v]);
    
    return b;
}

void texturePage::makeGLTexture(mesh* m, PtexTexture* r, GLenum format, GLenum typegl)
{
    if (faceIDs.size()<1) return;
    
    int bpp = r->numChannels()*Ptex::DataSize(r->dataType());
    
    GLuint fw = m->getFace(faceIDs[0])->ures;
    GLuint fh = m->getFace(faceIDs[0])->vres;
    
    char* pim = new char[fw*fh*faceIDs.size()*bpp];
    char* pimptr = pim;
    
    for (unsigned int i=0; i<faceIDs.size(); i++)
    {
        int fi = faceIDs[i];
        polygon* q = m->getFace(fi);
        if (q)
        {
            if (r->numFaces()==1 && r->meshType()==Ptex::mt_quad)
                r->getData(fi, pimptr+fw*(fh-1)*bpp, -fw*bpp);
            else 
                r->getData(fi, pimptr, 0);
        }
        pimptr += fw*fh*bpp;
    }
    if (!glIsTexture(id)) glGenTextures(1,&id);
    
    glBindTexture(GL_TEXTURE_2D_ARRAY,id);
    glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D_ARRAY,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage3D(GL_TEXTURE_2D_ARRAY,0,format,fw,fh,faceIDs.size(),0,format,typegl,pim);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    
    delete [] pim;
}

struct userOptions
{
    bool bgColorLight, showGridLines, envMode;
    
    userOptions()
    {
        bgColorLight = false;
        showGridLines = true;
        envMode = false;
    }
};

class PtexViewer
{
public:
    PtexViewer();
    bool load(std::string s, int faceid = -1, bool flat = false);
    void render(int selx = -1, int sely = -1);
    void orbit(int x, int y);
    void pan(int x, int y);
    void zoom(int z);
    void frameView(bool reset = false);
    void nextFace(int step);
    void firstLastFace(bool first);
    void allFaces();
    void toggleFullScreen();
    void menuAction(int i);
    void selectFace(int x, int y);
    void setDataScale(float s);
    void setDataBias(float b);
    void removeSceneData();
    void deleteShader();
    userOptions& options() { return _options; }
    
private:
    mesh _geometry;
    bool _mode3d, _envCube, _enable3D, _quads;
    BBox3f _bounds;
    camera _cam, _mainCam;
    std::string _filename, _typeStr;
    int _displayFace, _numFaces;
    userOptions _options;
    GLhandleARB _shaderProgram, _vertShader, _fragShader;
    
    void goToFace(int i);
    bool loadFile(std::string filename, int face = -1);
    void getWindowSize(double& width, double& height);
    void computeBounds();
    void renderChars(float x, float y, void *font, const char *string);
    void makeShader();
};

PtexViewer::PtexViewer()
{
    _displayFace = -1;
    _enable3D = true;
    _quads = true;
    _mainCam.setDistance(0);
    _shaderProgram = _vertShader = _fragShader = 0;
}

void PtexViewer::removeSceneData()
{
    _geometry.deleteData();
}

bool PtexViewer::load(std::string s, int faceid, bool flat)
{
    bool success = true;
    
    _enable3D = !flat;
    success = loadFile(s, faceid);
    if (success && faceid>=0) _displayFace = faceid;
    
    _filename = s;
    
    frameView(true);
    
    return success;
}

void PtexViewer::computeBounds()
{
    _bounds = _geometry.getBounds();
}

void PtexViewer::setDataScale(float s)
{
    glPixelTransferf(GL_RED_SCALE,s);
    glPixelTransferf(GL_GREEN_SCALE,s);
    glPixelTransferf(GL_BLUE_SCALE,s);
}

void PtexViewer::setDataBias(float b)
{
    glPixelTransferf(GL_RED_BIAS,b);
    glPixelTransferf(GL_GREEN_BIAS,b);
    glPixelTransferf(GL_BLUE_BIAS,b);
}

bool PtexViewer::loadFile(std::string filename, int face)
{
    Ptex::String ptexError;
    PtexPtr<PtexTexture> r(PtexTexture::open(filename.c_str(), ptexError, true));
    
    _quads = true;
    
    if (!r)
    {
        fprintf(stderr,"%s\n",ptexError.c_str());
        return false;
    }
    
    if (face>=r->numFaces())
    {
        std::cerr << "Invalid face ID requested" << std::endl;
        return false;
    }
    
    removeSceneData();
    
    _numFaces = r->numFaces();
    _typeStr = Ptex::DataTypeName(r->dataType());
    _quads = r->meshType()==Ptex::mt_quad;
    _mode3d = false;
    _envCube = false;
    
    if (_enable3D && r->numFaces()==6 && face<0)
    {
        // Test for environment cube case
        
        _envCube = true;
        
        int adjfaces[6][4] = {
            { 3, 5, 2, 4 }, // px
            { 3, 4, 2, 5 }, // nx
            { 4, 0, 5, 1 }, // py
            { 5, 0, 4, 1 }, // ny
            { 3, 0, 2, 1 }, // pz
            { 3, 1, 2, 0 }  // nz
        };
        
        int fi = 0;
        while (fi<6 && _envCube)
        {
            const Ptex::FaceInfo& f = r->getFaceInfo(fi);
            
            for (int v=0; v<4; v++)
                if (f.adjfaces[v]!=adjfaces[fi][v]) _envCube = false;
            
            fi++;
        }
    }
    
    if (_envCube)
    {
        // Construct environment cube geometry
        
        float cubeVerts[8][3] = { 
            {-1,-1,-1}, {1,-1,-1}, {1,-1, 1}, {-1,-1, 1},  
            {-1, 1,-1}, {1, 1,-1}, {1, 1, 1}, {-1, 1, 1} };

        int cubeInds[6][4] = {
            {2,1,5,6}, {0,3,7,4},
            {7,6,5,4}, {0,1,2,3},
            {3,2,6,7}, {1,0,4,5} };
        
        for (int i=0; i<8; i++) _geometry.addVertex(cubeVerts[i]);
        
        for (int fi=0; fi<6; fi++)
        {
            const Ptex::FaceInfo& f = r->getFaceInfo(fi);

            polygon q;
            for (int v=0; v<4; v++) q.indices[v] = cubeInds[fi][v];
            q.ures = f.res.u();
            q.vres = f.res.v();
            _geometry.addFace(q);
        }
        
        _mode3d = true;
    }
    else
    {
        if (_enable3D) _mode3d = _geometry.loadFromMetaData(r);
        if (_geometry.empty()) _geometry.makeFlatGeom(r);

        for (int fi=0; fi<r->numFaces(); fi++)
        {
            polygon* q = _geometry.getFace(fi);
            if (!q) continue;
            const Ptex::FaceInfo& f = r->getFaceInfo(fi);
            q->ures = f.res.u();
            q->vres = f.res.v();
        }
    }
    _geometry.makePages(r);
    
    computeBounds();
    
    return true;
}

void PtexViewer::getWindowSize(double& width, double& height)
{
    GLint viewport[4];  
    glGetIntegerv(GL_VIEWPORT, viewport);
    width = ((double)viewport[2]);
    height = ((double)viewport[3]);
}

void PtexViewer::orbit(int x, int y)
{
    if (_mode3d && _displayFace<0) _cam.rotateXY(x,-y);
}

void PtexViewer::pan(int x, int y)
{
    double wsize[2];
    getWindowSize(wsize[0], wsize[1]);
    wsize[0] /= 2.0;
    wsize[1] /= 2.0;
    
    float winasp = wsize[0]/wsize[1];
  
    double panX, panY;
    double panscale = _cam.getDistance()*0.5;
    panX = panscale*((double)(x))/wsize[0];
    panY = panscale*((double)(-y))/(wsize[1]*winasp);
    
    _cam.pan(panX,panY,(_mode3d&&_displayFace<0));
}

void PtexViewer::zoom(int z)
{
    _cam.zoom(z);
}

void PtexViewer::frameView(bool reset)
{
    if (reset) _cam.reset();
    
    if (_displayFace>=0)
    {
        _cam.setLookAt(Vec3f(0.5, 0.5, 0));
        _cam.setDistance(1.41);
    }
    else
    {
        if (_mode3d && _envCube && reset)
        {
            _cam.setLookAt(Vec3f(0,0,0));
            _cam.setDistance(1e-4);
        }
        else
        {
            _cam.setLookAt(_bounds.getCenter());
            _cam.setDistance(_bounds.diagonal());
        }
    }
}

void PtexViewer::renderChars(float x, float y, void* font, const char* string)
{
    glRasterPos3f(x, y, 0);
    for (const char* c=string; *c != '\0'; c++) glutBitmapCharacter(font, *c);
}

void PtexViewer::selectFace(int x, int y)
{
    if (_numFaces<2 || _displayFace>=0) return;
    
    GLuint* selbuf = new GLuint[_numFaces*4];
    glSelectBuffer(_numFaces*4,selbuf);
    
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);
    render(x,y);
    glPopName();
    GLint numhits = glRenderMode(GL_RENDER);

    GLuint *ptr,names,nameID;
    float z1;
    GLdouble dist = 1.0;
    ptr = (GLuint*)selbuf;
    int foundFace = -1;
    
    for (int i=0; i<numhits; i++)
    {
        names = *ptr;
        ptr++;
        z1 = (float) *ptr / 0xffffffff;
        ptr++;
        ptr++;
        nameID = *ptr;
        ptr++;
        for (unsigned int n=1; n<names; n++) ptr++;
        if (nameID == INVALIDNAME ) continue;
        if (z1<=dist)
        {
            dist = z1;
            foundFace = nameID;
        }
    }
    
    delete [] selbuf;
    
    if (foundFace>=0 && foundFace<_numFaces)
    {
        if (_displayFace<0) _mainCam = _cam;
        
        _displayFace = foundFace;

        frameView();
    }
}
        
void PtexViewer::render(int selx, int sely)
{
    if (!_shaderProgram) makeShader();

    if (_options.bgColorLight) glClearColor(0.8,0.8,0.8,1);
    else glClearColor(0.0,0.0,0.0,1);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    
    GLint vp[4];  
    glGetIntegerv(GL_VIEWPORT, vp);
    float width = vp[2];
    float height = vp[3];
    float winasp = width/height;
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    bool selectMode = selx>=0 && sely>=0;
    
    if (selectMode)
        gluPickMatrix((double)selx,(double)vp[3]-sely,3,3,vp);
    
    if (!_mode3d||_displayFace>=0)
    {
        float sc = _cam.getDistance()/2;
        glOrtho(-sc, sc, -sc/winasp, sc/winasp, -100, 100);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(-_cam.getLookAt()[0],-_cam.getLookAt()[1],-_cam.getLookAt()[2]);
    }
    else
    {
        _cam.applyProjectionTransform(winasp, _bounds);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        _cam.applyModelViewTransform();
    }
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0,1.0);
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glUseProgramObjectARB(_shaderProgram);
    _geometry.draw(_displayFace);
    glUseProgramObjectARB(0);
    if (!selectMode && _options.showGridLines)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        glColor4f(0.0, 0.0, 0.0, 0.25);
        _geometry.draw(_displayFace);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDisable(GL_BLEND);
        glDisable(GL_LINE_SMOOTH);
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_DEPTH_TEST);
    
    if (selectMode) return;
    
    // Display text information
    char str[64];
    
    if (_options.bgColorLight) glColor3f(0.0, 0.0, 0.0);
    else glColor3f(0.9, 0.9, 0.9);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    if (!_options.envMode)
    {
        if (_numFaces==1) strcpy(str, "1 face");
        else sprintf(str,"%i faces",_numFaces);

        renderChars(-0.8,-0.92,GLUT_BITMAP_8_BY_13,str);
        if (_displayFace>=0 && _numFaces!=1)
        {
            sprintf(str,"Face ID: %i",_displayFace);
            renderChars(-0.8,0.92,GLUT_BITMAP_8_BY_13,str);
        }
        if (_numFaces==1)
        {
            int ures = _geometry.getFace(0)->ures;
            int vres = _geometry.getFace(0)->vres;
            sprintf(str,"Res: %i x %i",ures,vres);
            renderChars(-0.8,0.92,GLUT_BITMAP_8_BY_13,str);
        }
        else if (_displayFace>=0)
        {
            int ures = _geometry.getFace(_displayFace)->ures;
            int vres = _geometry.getFace(_displayFace)->vres;
            sprintf(str,"Face Res: %i x %i",ures,vres);
            renderChars(-0.2,0.92,GLUT_BITMAP_8_BY_13,str);
        }
        renderChars(0.25,-0.92,GLUT_BITMAP_8_BY_13,(char*)_typeStr.c_str());
    }
}

void PtexViewer::goToFace(int i)
{
    if (i>=0 && _displayFace<0) _mainCam = _cam;
    else if (i<0) _cam = _mainCam;
    
    _displayFace = i;
    
    if (_displayFace>=0 || (i<0&&_mainCam.getDistance()==0)) frameView();
}

void PtexViewer::nextFace(int step)
{
    if (_numFaces==1) return;
    if (step<0 && _displayFace<0) return;
    
    goToFace(std::min(std::max(_displayFace+step,-1),_numFaces-1));
}

void PtexViewer::firstLastFace(bool first)
{
    if (_numFaces==1) return;
    
    goToFace(first ? 0 : _numFaces-1);
}

void PtexViewer::allFaces()
{
    if (_numFaces==1 || _displayFace<0) return;
    
    goToFace(-1);
}

void PtexViewer::toggleFullScreen()
{
    static bool fullscreen = false;
    static int win[2];
    if (!fullscreen)
    {
        win[0] = glutGet(GLUT_WINDOW_WIDTH);
        win[1] = glutGet(GLUT_WINDOW_HEIGHT);
        glutFullScreen();
        fullscreen = true;
    }
    else
    {
        glutReshapeWindow(win[0],win[1]);
        fullscreen = false;
    }
}

void PtexViewer::menuAction(int i)
{
    switch (i)
    {
        case MENU_BGCOLOR :
            _options.bgColorLight = !_options.bgColorLight;
            break;
        case MENU_GRIDLINES :
            _options.showGridLines = !_options.showGridLines;
            break;
        case MENU_RESETVIEW :
            frameView(true);
            break;
        case MENU_FOCUSVIEW :
            frameView();
            break;
        case MENU_FULLSCREEN :
            toggleFullScreen();
            break;
        case MENU_QUIT :
            exit(EXIT_SUCCESS);
            break;
    }
}

void PtexViewer::makeShader()
{
    deleteShader();
    
    const char *vs =
    "void main( void ) { gl_Position = ftransform(); gl_TexCoord[0] = gl_MultiTexCoord0; }";

    const char *fs = 
    "varying vec4 gl_TexCoord[4]; \
    uniform sampler2DARRAY tex; \
    void main( void ) { \
    gl_FragColor = tex2DARRAY(tex,gl_TexCoord[0].xyz);\
    }";

    _vertShader  = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
    _fragShader  = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
    
    glShaderSourceARB(_vertShader, 1, &vs,NULL);
    glShaderSourceARB(_fragShader, 1, &fs,NULL);

    glCompileShaderARB(_vertShader);
    glCompileShaderARB(_fragShader);

    _shaderProgram = glCreateProgramObjectARB();
    glAttachObjectARB(_shaderProgram,_vertShader);
    glAttachObjectARB(_shaderProgram,_fragShader);

    glLinkProgramARB(_shaderProgram);
    
    glUseProgramObjectARB(0);
}

void PtexViewer::deleteShader()
{
    if (_vertShader)
    {
        if (_shaderProgram) glDetachObjectARB(_shaderProgram,_vertShader);
        glDeleteObjectARB(_vertShader);
    }
    if (_fragShader)
    {
        if (_shaderProgram) glDetachObjectARB(_shaderProgram,_fragShader);
        glDeleteObjectARB(_fragShader);
    }
    if (_shaderProgram) glDeleteObjectARB(_shaderProgram);
}

struct inputState
{
    int modifier;
    int mbase[2];
    bool buttons[5];
};

PtexViewer w;
inputState input;

void renderScene()
{
    w.render();
    glutSwapBuffers();
}

void closeEvent()
{
    w.removeSceneData();
    w.deleteShader();
}

void resizeEvent(int width, int height)
{
    glViewport(0,0,width,height);
}

void mouseButtonEvent(int button, int state, int x, int y)
{
    input.modifier = glutGetModifiers();
    if (button<5) input.buttons[button] = !state;
    
    if (button==3)
    {
        w.zoom(10);
        glutPostRedisplay();
    }
    else if (button==4)
    {
        w.zoom(-10);
        glutPostRedisplay();
    }
    
    if (w.options().envMode) return;
    
    if (input.modifier == 0 && state==0)
    {
        if (input.buttons[0]&&!input.buttons[1]&&!input.buttons[2])
        {
            w.selectFace(x,y);
            glutPostRedisplay();
        }
    }
}

void mouseEvent(int x, int y)
{
    if (input.modifier == GLUT_ACTIVE_ALT || w.options().envMode)
    {
        if (input.buttons[0]&&!input.buttons[1]&&!input.buttons[2])
        {
            w.orbit( x-input.mbase[0], y-input.mbase[1] );
        }
        else if (!input.buttons[0]&&input.buttons[1]&&!input.buttons[2])
        {
            w.pan( x-input.mbase[0], y-input.mbase[1] );
        }
        else if ((!input.buttons[0]&&!input.buttons[1]&&input.buttons[2])||
                 (input.buttons[0]&&input.buttons[1]&&!input.buttons[2]))
        {
            w.zoom( x-input.mbase[0] );
        }
    }
    
    input.mbase[0] = x;
    input.mbase[1] = y;
    glutPostRedisplay();
}

void mousePassiveEvent(int x, int y)
{
    input.mbase[0] = x;
    input.mbase[1] = y;
}

void keyboardEvent(unsigned char key, int /*x*/, int /*y*/)
{
    if (key=='f') w.frameView();
    if (key=='r') w.frameView(true);
    if (key==27)
    {
        if (w.options().envMode) w.toggleFullScreen();
        else w.allFaces();
    }
    if (key=='q') exit(EXIT_SUCCESS);
    glutPostRedisplay();
}

void keySpecEvent(int key, int /*x*/, int /*y*/)
{
    switch (key)
    {
        case GLUT_KEY_UP: w.allFaces(); break;
        case GLUT_KEY_RIGHT: w.nextFace(1); break;
        case GLUT_KEY_LEFT: w.nextFace(-1); break;
        case GLUT_KEY_PAGE_UP: w.nextFace(100); break;
        case GLUT_KEY_PAGE_DOWN: w.nextFace(-100); break;
        case GLUT_KEY_HOME: w.firstLastFace(true); break;
        case GLUT_KEY_END: w.firstLastFace(false); break;
    }
    glutPostRedisplay();
}

void processMenuEvents(int option)
{
    w.menuAction(option);
    glutPostRedisplay();
}

void createMenus()
{
    glutCreateMenu(processMenuEvents);
    glutAddMenuEntry("Toggle Background",MENU_BGCOLOR);
    glutAddMenuEntry("Toggle Edge Lines",MENU_GRIDLINES);
    glutAddMenuEntry("Reset View",MENU_RESETVIEW);
    glutAddMenuEntry("Focus View",MENU_FOCUSVIEW);
    glutAddMenuEntry("Full Screen",MENU_FULLSCREEN);
    glutAddMenuEntry("Quit",MENU_QUIT);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void showHelp(bool full = false)
{
    std::cout << "Usage: ptxview [options] filename.ptx\n";
    
    if (full)
    {
        std::cout << "\n";
        std::cout << "Options:\n";
        std::cout << "        -help             print this message\n";
        std::cout << "        -face n           display given face ID\n";
        std::cout << "        -flat             display face data flat (2D viewing)\n";
        std::cout << "        -scale s          scale image data by s\n";
        std::cout << "        -bias b           bias image data by b\n";
        std::cout << "\n";
        std::cout << "Navigation:\n";
        std::cout << "        Alt-Button 1:     rotate scene (3D mode only)\n";
        std::cout << "        Alt-Button 2:     pan scene\n";
        std::cout << "        Alt-Buttons 1&2:  zoom scene\n";
        std::cout << "\n";
        std::cout << "Hot Keys:\n";
        std::cout << "        Right Arrow:      view next face\n";
        std::cout << "        Left Arrow:       view previous face\n";
        std::cout << "        Page Up:          advance 100 faces\n";
        std::cout << "        Page Down:        go back 100 faces\n";
        std::cout << "        Home:             view first face\n";
        std::cout << "        Page Down:        view last face\n";
        std::cout << "        Up Arrow:         view all faces\n";
        std::cout << "        Escape:           view all faces\n";
        std::cout << "        r:                reset view\n";
        std::cout << "        f:                center data on screen\n";
        std::cout << "        q:                quit\n";
    }
    
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc<2) showHelp();
    
    std::string filename = "";
    int faceid = -1;
    bool flat = false, fs = false;
    float datascale = 0, databias = 0;
    
    int i = 1;
    while (i<argc)
    {
        std::string str(argv[i]);
        if (str=="-face") faceid = atoi(argv[++i]);
        else if (str=="-flat") flat = true;
        else if (str=="-scale")
        {
            float f = atof(argv[++i]);
            if (!isnan(f) && f!=0) datascale = f;
        }
        else if (str=="-bias")
        {
            float f = atof(argv[++i]);
            if (!isnan(f)) databias = f;
        }
        else if (str=="-fullscreen")
        {
            fs = true;
        }
        else if (str=="-help" || str=="-h") showHelp(true);
        else filename = str;
        i++;
    }
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(128, 128);
    glutInitWindowSize(512, 512);
    std::string s = filename;
    if (s.length()>=64)
    {
        size_t pos = s.rfind('/');
        if (pos!=std::string::npos && s.length()-pos<=64)
            s = s.substr(pos+1,s.length()-pos-1);
        else
            s = s.substr(s.length()-64,64);
    }

    std::string windowtitle = "ptxview ";
    windowtitle += s;
    if (windowtitle.size() > 128)
        windowtitle.resize(128);

    glutCreateWindow(windowtitle.c_str());

    glutDisplayFunc(renderScene);
    glutCloseFunc(closeEvent);
    glutReshapeFunc(resizeEvent);
    glutMotionFunc(mouseEvent);
    glutMouseFunc(mouseButtonEvent);
    glutPassiveMotionFunc(mousePassiveEvent);
    glutKeyboardFunc(keyboardEvent);
    glutSpecialFunc(keySpecEvent);

    glewInit();

    createMenus();
    
    if (datascale!=0) w.setDataScale(datascale);
    if (databias!=0) w.setDataBias(databias);
    
    if (filename=="") showHelp();
    if (!w.load(filename, faceid, flat)) exit(1);
    
    if (filename.find(".penv")!=std::string::npos && !flat && faceid<0)
    {
        w.options().envMode = true;
        w.options().showGridLines = false;
    }
    
    if (fs) w.toggleFullScreen();
    
    glutMainLoop();
    
    return 0;
}
