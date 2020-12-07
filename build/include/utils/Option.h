#ifndef OPTION_H
#define OPTION_H

namespace Utils {

    template<typename T>
        class Option {
            public:
                virtual const bool has_value() const = 0;
                virtual ~Option() = 0;
        };
    template<typename T>
    Option<T>::~Option() {}

    template<typename T>
        const Option<T> &None();

    template<typename T>
        class Some : public Option<T> {
            public:
                Some(const T &v) : val(v) {}
                virtual ~Some() {}
                virtual const bool has_value() const { return true; };
                const T val;
        };

    template<typename T>
        class NoneType : public Option<T> {
            template<typename U> 
                friend const Option<U> &None();
            public:
            virtual ~NoneType() {}
            virtual const bool has_value() const { return false; };
            private:
            NoneType() {}
        };

    template<typename T>
        const Option<T> &None() {
            static NoneType<T> x;
            return x;
        }

    template<typename T>
        bool operator==(const Option<T> &lhs, const Option<T> &rhs) {
            if (lhs.has_value() != rhs.has_value()) { // one is Some, other is None
                return false;
            } else if (!lhs.has_value()) { // both are None
                return true;
            } else {
                // Now it should be safe
                return ((static_cast<const Some<T> &>(lhs)).val == (static_cast<const Some<T> &>(rhs)).val);
            }
        }

    template<typename T>
        bool operator!=(const Option<T> &lhs, const Option<T> &rhs) {
            return !(lhs == rhs);
        }

}
#endif
