//Zachariah Vicars
//Provides some useful typedefs
#pragma once
#include <vector>
#include <array>
#include <string>
#include <map>

template <class T> 
using Vec = std::vector<T>;

template <class T> 
using Vec2 = std::array<T,2>;

template <class T> 
using Vec3 = std::array<T,3>;

template <class T> 
using Vec4 = std::array<T,4>;

#ifndef USE_FLOAT
using real = double;
#else
using real = float;
#endif

using string = std::string;

template <class T>
using strmap = std::map<std::string, T>;

using size_t = std::size_t;
