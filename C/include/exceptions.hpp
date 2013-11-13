#ifndef turb_exceptions_h
#define turb_exceptions_h

#include <string>
#include <sstream>

namespace turb {

    struct Exception {
    protected:
        std::string message;

    public:
        Exception() { message = ""; }
        Exception(std::string msg) { message = msg; }
        Exception(std::stringstream *output) { message = output->str();}

        std::string what() const { return message;}
    };

    struct InstanceNotFound : public Exception {
        using Exception::Exception;
    };

    struct ProgramDeathRequest : public Exception {
        using Exception::Exception;
    };

    struct RemoveOutput : public Exception {
        using Exception::Exception;
    };
}

#endif
