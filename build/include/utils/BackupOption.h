#ifndef OPTION_H
#define OPTION_H

#warning "Library code, remove this include when finished"
#include<iostream>

namespace Utils {

    template<typename T>
        class Option {
            public:
                virtual ~Option() = 0;
        };
    template<typename T>
    Option<T>::~Option() {}

    template<typename T>
        const Option<T> &None();

    template<typename T>
        const Option<T> &Some(const T &val);

    template<typename T>
        class SomeType : public Option<T> {
            template<typename U>
                friend const Option<U> &Some(const U &val);
            public:
            virtual ~SomeType() {}
            const T val;
            private:
            SomeType(const T &v) : val(v) {}
        };

    template<typename T>
        class NoneType : public Option<T> {
            template<typename U> 
                friend const Option<U> &None();
            public:
            virtual ~NoneType() {}
            private:
            NoneType() {}
        };

    template<typename T>
        const Option<T> &None() {
            static NoneType<T> x;
            return x;
        }

    template<typename T>
        const Option<T> &Some(const T &val) {
            return SomeType<T>(val);
        }

    template<typename T>
        bool operator==(NoneType<T> &lhs, Option<T> &rhs) {
            return false;
        }
    template<typename T>
        bool operator==(Option<T> &lhs, NoneType<T> &rhs) {
            return false;
        }
    template<typename T>
        bool operator==(NoneType<T> &lhs, NoneType<T> &rhs) {
            return true;
        }
    template<typename T>
        bool operator==(const Option<T> &lhs, const Option<T> &rhs) {
            // Now it should be safe
            std::cout << "Correct operator is called" << std::endl;
            std::cout << (static_cast<const SomeType<T> &>(lhs)).val << " " << (static_cast<const SomeType<T> &>(rhs)).val << std::endl;
            return ((static_cast<const SomeType<T> &>(lhs)).val == (static_cast<const SomeType<T> &>(rhs)).val);
        }
    template<typename T>
        bool operator!=(const Option<T> &lhs, const Option<T> &rhs) {
            return !(lhs == rhs);
        }

}
#endif
