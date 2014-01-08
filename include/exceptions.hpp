#ifndef turb_exceptions_h
#define turb_exceptions_h

#include <string>
#include <sstream>

namespace turb {

    /*! \brief A basis of all turb exceptions.
     *
     *  Most use cases require just a simple
     *  inheritance.
     */
    struct Exception {
    protected:
        /// Info message the exception provides
        std::string message;

    public:
        Exception() { message = ""; }
        /*! \brief Sets the info message to something
         *      meaningful.
         *
         *  A std::stringstream accepting form is also
         *  provided.
         */
        Exception(std::string msg) { message = msg; }
        Exception(std::stringstream *output) { message = output->str();}
        Exception(double num) { std::stringstream tmp; tmp << num; message = tmp.str();}

        /*! \brief Returns the info message
         *  \return info message
         */
        std::string what() const { return message;}
    };

    /*! \brief Called when an instance is requested that
     *      could not be found at the moment.
     */
    struct InstanceNotFound : public Exception {
        //using Exception::Exception;
    public:
        InstanceNotFound() { message = ""; }
        InstanceNotFound(std::string msg) { message = msg; }
        InstanceNotFound(std::stringstream *output) { message = output->str();}
        InstanceNotFound(double num) { std::stringstream tmp; tmp << num; message = tmp.str();}
    };

    /*! \brief Indicates a fatal error requiring
     *      immediate program termination.
     */
    struct ProgramDeathRequest : public Exception {
        //using Exception::Exception;
    public:
        ProgramDeathRequest() { message = ""; }
        ProgramDeathRequest(std::string msg) { message = msg; }
        ProgramDeathRequest(std::stringstream *output) { message = output->str();}
        ProgramDeathRequest(double num) { std::stringstream tmp; tmp << num; message = tmp.str();}
    };

    /*! \brief Thrown when currently used output file
     *      needs to be deleted properly.
     */
    struct RemoveOutput : public Exception {
        //using Exception::Exception;
    public:
        RemoveOutput() { message = ""; }
        RemoveOutput(std::string msg) { message = msg; }
        RemoveOutput(std::stringstream *output) { message = output->str();}
        RemoveOutput(double num) { std::stringstream tmp; tmp << num; message = tmp.str();}
    };

    /*! \brief Thrown if there is no result available
     *      to be thrown.
     */
    struct NoResult : public Exception {
        //using Exception::Exception;
    public:
        NoResult() { message = ""; }
        NoResult(std::string msg) { message = msg; }
        NoResult(std::stringstream *output) { message = output->str();}
        NoResult(double num) { std::stringstream tmp; tmp << num; message = tmp.str();}
    };

    struct OutOfBounds : public Exception {
        //using Exception::Exception;
    public:
        OutOfBounds() { message = ""; }
        OutOfBounds(std::string msg) { message = msg; }
        OutOfBounds(std::stringstream *output) { message = output->str();}
        OutOfBounds(double num) { std::stringstream tmp; tmp << num; message = tmp.str();}
    };
}

#endif
