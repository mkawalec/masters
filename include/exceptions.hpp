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
        using Exception::Exception;
    };

    /*! \brief Indicates a fatal error requiring 
     *      immediate program termination.
     */
    struct ProgramDeathRequest : public Exception {
        using Exception::Exception;
    };

    /*! \brief Thrown when currently used output file
     *      needs to be deleted properly.
     */
    struct RemoveOutput : public Exception {
        using Exception::Exception;
    };

    /*! \brief Thrown if there is no result available
     *      to be thrown.
     */
    struct NoResult : public Exception {
        using Exception::Exception;
    };

    struct OutOfBounds : public Exception {
        using Exception::Exception;
    };
}

#endif
