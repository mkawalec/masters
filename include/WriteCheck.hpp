#ifndef turb_WriteCheck_h
#define turb_WriteCheck_h

namespace turb {

    template <typename T1, typename T2>
    class WriteCheck {
    private:
        T1* store;
        T2* value;
        int index;

    public:
        WriteCheck(T1* s, T2* v, int i) : store(s), value(v), index(i) {}

        WriteCheck& operator=(T2 const& rhs)
        {
            *value = rhs;
            store->update_state(index);
            return *this;
        }

        WriteCheck& operator-=(T2 const& rhs)
        {
            *value -= rhs;
            store->update_state(index);
            return *this;
        }

        WriteCheck& operator+=(T2 const& rhs)
        {
            *value += rhs;
            store->update_state(index);
            return *this;
        }

        operator T2 const&()
        {
            return *value;
        }
    };
}

#endif

