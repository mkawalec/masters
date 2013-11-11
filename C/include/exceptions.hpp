#ifndef turb_exceptions_h
#define turb_exceptions_h

#include <string>

namespace turb {

    struct Exception {
    protected:
        std::string message;

    public:
        std::string what() { return message;}
        Exception(std::string msg) { message = msg; }
        Exception() { message = ""; }
    };

    struct SerializerNotFound : public Exception {
        using Exception::Exception;
    };
}

#endif
