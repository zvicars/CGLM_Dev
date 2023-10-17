#pragma once
#include <fstream>
#include <iostream>
#include "typedefs.hpp"
#include "Matrix.hpp"
#include "Assert.hpp"
static void binary_bool_write(std::ofstream& fout, const std::vector<bool>& x, Vec3<std::size_t>& size)
{
    std::size_t n = size[0]*size[1]*size[2];
    std::size_t sx = size[0], sy = size[1], sz = size[2];
    fout.write((char*)&sx, sizeof(std::size_t));
    fout.write((char*)&sy, sizeof(std::size_t));
    fout.write((char*)&sz, sizeof(std::size_t));
    for(std::vector<bool>::size_type i = 0; i < n;)
    {
        unsigned char aggr = 0;
        for(unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1)
            if(x.at(i))
                aggr |= mask;
        fout.write((const char*)&aggr, sizeof(unsigned char));
    }
    return;
}

static void binary_bool_read(std::ifstream& fin, std::vector<bool>& x, Vec3<std::size_t>& size)
{
    std::size_t sx, sy, sz;
    fin.read((char*)&sx, sizeof(std::size_t));
    if(fin.fail()) return;
    fin.read((char*)&sy, sizeof(std::size_t));
    if(fin.fail()) return;
    fin.read((char*)&sz, sizeof(std::size_t));
    if(fin.fail()) return;
    size[0] = sx;
    size[1] = sy;
    size[2] = sz;
    x.resize(size[0]*size[1]*size[2]);
    for(std::vector<bool>::size_type i = 0; i < x.size();)
    {
        unsigned char aggr;
        fin.read((char*)&aggr, sizeof(unsigned char));
        if(fin.fail()) return;
        for(unsigned char mask = 1; mask > 0 && i < x.size(); ++i, mask <<= 1)
            x.at(i) = aggr & mask;
    }
    return;
}

static void binary_char_write(std::ofstream& fout, const std::vector<char>& x, Vec3<std::size_t>& size)
{
    std::size_t n = size[0]*size[1]*size[2];
    std::size_t sx = size[0], sy = size[1], sz = size[2];
    fout.write((char*)&sx, sizeof(std::size_t));
    fout.write((char*)&sy, sizeof(std::size_t));
    fout.write((char*)&sz, sizeof(std::size_t));
    for(std::vector<char>::size_type i = 0; i < n; i++)
    {
        char x_i = x[i];
        fout.write((const char*)&x_i, sizeof(char));
    }
    return;
}

static void binary_char_read(std::ifstream& fin, std::vector<char>& x, Vec3<std::size_t>& size)
{
    std::size_t sx, sy, sz;
    fin.read((char*)&sx, sizeof(std::size_t));
    if(fin.fail()) return;
    fin.read((char*)&sy, sizeof(std::size_t));
    if(fin.fail()) return;
    fin.read((char*)&sz, sizeof(std::size_t));
    if(fin.fail()) return;
    size[0] = sx;
    size[1] = sy;
    size[2] = sz;
    x.resize(size[0]*size[1]*size[2]);
    for(std::vector<char>::size_type i = 0; i < x.size(); i++)
    {
        fin.read((char*)&x[i], sizeof(char));
        if(fin.fail()) return;
    }
    return;
}

static void binary_real_write(std::ofstream& fout, const std::vector<real>& x, Vec3<std::size_t>& size)
{
    std::size_t n = size[0]*size[1]*size[2];
    std::size_t sx = size[0], sy = size[1], sz = size[2];
    fout.write((char*)&sx, sizeof(std::size_t));
    fout.write((char*)&sy, sizeof(std::size_t));
    fout.write((char*)&sz, sizeof(std::size_t));
    for(std::vector<real>::size_type i = 0; i < n; i++)
    {
        real x_i = x[i];
        fout.write((const char*)&x_i, sizeof(real));
    }
    return;
}

static void binary_real_read(std::ifstream& fin, std::vector<real>& x, Vec3<std::size_t>& size)
{
    std::size_t sx, sy, sz;
    fin.read((char*)&sx, sizeof(std::size_t));
    if(fin.fail()) return;
    fin.read((char*)&sy, sizeof(std::size_t));
    if(fin.fail()) return;
    fin.read((char*)&sz, sizeof(std::size_t));
    if(fin.fail()) return;
    size[0] = sx;
    size[1] = sy;
    size[2] = sz;
    x.resize(size[0]*size[1]*size[2]);
    for(std::vector<real>::size_type i = 0; i < x.size(); i++)
    {
        fin.read((char*)&x[i], sizeof(real));
        if(fin.fail()) return;
    }
    return;
}
//this one is used for reading reals, must save/load with the same precision/endianness
static bool readFileIntoArray(std::ifstream& fin, const Vec3<std::size_t>& size, Matrix3d<real>& data){
    std::vector<real> f_in;
    Vec3<std::size_t> size_read;
    std::vector<real> data_read;
    binary_real_read(fin, data_read, size_read);
    FANCY_ASSERT(size_read[0] == size[0] && size_read[1] == size[1] && size_read[2] == size[2], "sizes do not match");
    data.initialize(size_read);
    for(std::size_t i = 0; i < data_read.size(); i++){
        Vec3<std::size_t> pos = data.map1N(i);
        auto new_idx = data.mapN1(pos);
        data.at(pos) = data_read[i];
        //std::cout << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << data.at(pos) << std::endl;
        //std::cin.get();
    }
    return 0;
}
//this one is used for the array of states, uses char array internally because bool arrays can behave oddly
static bool readFileIntoArray(std::ifstream& fin, const Vec3<std::size_t>& size, Matrix3d<char>& data){
    std::vector<real> f_in;
    Vec3<std::size_t> size_read;
    std::vector<bool> data_read;
    binary_bool_read(fin, data_read, size_read);
    FANCY_ASSERT(size_read[0] == size[0] && size_read[1] == size[1] && size_read[2] == size[2], "sizes do not match");
    data.initialize(size_read);
    for(int i = 0; i < data_read.size(); i++){
        data.at_1d(i) = data_read[i];
    }
    return 0;
}