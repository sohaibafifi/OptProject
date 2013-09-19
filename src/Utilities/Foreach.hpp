// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France
// This file is influenced by the Qt foreach macro and is not completely tested yet


// FIXME : test this macro
#pragma once
namespace Utilities{

#if defined(CC_GNU) && !defined(CC_INTEL) && !defined(CC_RVCT)
template <typename T>
class ForeachContainer {
public:
    inline ForeachContainer(const T& t) : c(t), brk(0), i(c.begin()), e(c.end()) { }
    const T c;
    int brk;
    typename T::const_iterator i, e;
};

#define FOREACH(variable, container)                                \
for (ForeachContainer<__typeof__(container)> _container_(container); \
     !_container_.brk && _container_.i != _container_.e;              \
     __extension__  ({ ++_container_.brk; ++_container_.i; }))                       \
    for (variable = *_container_.i;; __extension__ ({--_container_.brk; break;}))
#else

struct ForeachContainerBase {};

template <typename T>
class ForeachContainer : public ForeachContainerBase {
public:
    inline ForeachContainer(const T& t): c(t), brk(0), i(c.begin()), e(c.end()){}
    const T c;
    mutable int brk;
    mutable typename T::const_iterator i, e;
    inline bool condition() const { return (!brk++ && i != e); }
};

template <typename T> inline T *_ForeachPointer(const T &) { return 0; }

template <typename T> inline ForeachContainer<T> _ForeachContainerNew(const T& t)
{ return ForeachContainer<T>(t); }

template <typename T>
inline const ForeachContainer<T> *_ForeachContainer(const ForeachContainerBase *base, const T *)
{ return static_cast<const ForeachContainer<T> *>(base); }
#if defined(CC_DIAB)
// VxWorks DIAB generates unresolvable symbols, if container is a function call
#  define FOREACH(variable,container)                                                             \
    if(0){}else                                                                                     \
    for (const ForeachContainerBase &_container_ = _ForeachContainerNew(container);                \
         _ForeachContainer(&_container_, (__typeof__(container) *) 0)->condition();       \
         ++_ForeachContainer(&_container_, (__typeof__(container) *) 0)->i)               \
        for (variable = *_ForeachContainer(&_container_, (__typeof__(container) *) 0)->i; \
             _ForeachContainer(&_container_, (__typeof__(container) *) 0)->brk;           \
             --_ForeachContainer(&_container_, (__typeof__(container) *) 0)->brk)

#else
#  define FOREACH(variable, container) \
    for (const ForeachContainerBase &_container_ = _ForeachContainerNew(container); \
         _ForeachContainer(&_container_, true ? 0 : _ForeachPointer(container))->condition();       \
         ++_ForeachContainer(&_container_, true ? 0 : _ForeachPointer(container))->i)               \
        for (variable = *_ForeachContainer(&_container_, true ? 0 : _ForeachPointer(container))->i; \
             _ForeachContainer(&_container_, true ? 0 : _ForeachPointer(container))->brk;           \
             --_ForeachContainer(&_container_, true ? 0 : _ForeachPointer(container))->brk)
#endif // MSVC6 || MIPSpro

#endif

#define FOREVER for(;;)
#ifndef foreach
#  define foreach FOREACH
#endif
#ifndef forever
#  define forever FOREVER
#endif
}// end namespace
