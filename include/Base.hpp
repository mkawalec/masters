#ifndef turb_Base_h
#define turb_Base_h

#include "helpers.hpp"

#include <list>
#include <string>

namespace turb {

    /*! \brief Provides a common scafold for dynamic worker
     *      classes.
     *
     *  Inherit from this class if you wish to add a new
     *  worker type. Functions for easy usage of dynamically
     *  allocated worker instances are provided
     */
    template <typename T>
    class Base {
    protected:
        /*! \brief Called when an instance wishes to stop
         *      being visible to possible users.
         */
        void unregister();

        /// Contains currently available worker classes
        static std::list<T*> available;

        /*! \brief Gives additional information about a worker.
         *  \return the additional information
         *
         *  Usually used to provide description of
         *  class-specific command line parameters
         */
        virtual std::string additional_info() { return "";}

    public:
        /*! \brief Called to request a worker class
         *      with a specified name.
         *
         *  \param name name of a requested worker class
         *  \exception InstanceNotFound if no worker with a
         *      specified name is found
         *
         *  \return pointer to instance OR NULL when name is
         *      an empty std::string
         */
        static T* choose(std::string name);

        /*! \brief Provides a human-readable information
         *      about currently available workers.
         *
         *  \return available workers description
         */
        static std::string list_available();

        /*! \brief Parses class-specific command line parameters.
         *  \param argc number of command-line parameters
         *  \param argv array of command-line parameters
         */
        virtual void parse_params(int argc, const char *argv[])
            { unused(argc); unused(argv);}

        /// Name of the current class, specifying its function
        std::string name;

        /// Human-readable action description
        std::string description;

        /*! Name of the type of class the particular class
         *  represents (for instance "computer" etc.)
         */
        static std::string class_name;

        /*! \brief Used to obtain clones for polymorphic
         *      classes inheriting from Base<T>
         *
         *  \return pointer to new cloned instance
         */
        virtual Base<T>* clone() const = 0;
    };
}

#include "Base.cpp"

#endif
