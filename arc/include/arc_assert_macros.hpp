/*********************************************************************************
 *  ARC-SRTK - Single Frequency RTK Pisitioning Library
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Created on: July 08, 2017
 *      Author: SuJingLan
 *********************************************************************************/

/**
 * @file arc_assert_macros.hpp
 * @brief This file contains some useful assert macros.
 * @author SuJingLan
 */

#ifndef ARC_ASSERT_MACROS_HPP
#define ARC_ASSERT_MACROS_HPP

#include <stdexcept>
#include <sstream>
#include <typeinfo>
#include "arc_source_file_pos.hpp"

//! Macro for defining an exception with a given parent
/// (std::runtime_error should be top parent)
/// adapted from ros/drivers/laser/hokuyo_driver/hokuyo.h
#define ARC_DEFINE_EXCEPTION(exceptionName, exceptionParent)				\
  class exceptionName : public exceptionParent {						    \
public: 																    \
  exceptionName(const char * message) : exceptionParent(message) {}		    \
  exceptionName(std::string const & message) : exceptionParent(message) {}  \
  virtual ~exceptionName() throw() {}									    \
};

/// \brief ARC Main namespace of this package.
namespace ARC {

    namespace detail {

        template<typename ARC_EXCEPTION_T>
        inline void ARC_throw_exception(std::string const & exceptionType, ARC::source_file_pos sfp, std::string const & message)
        {
            std::stringstream ARC_assert_stringstream;
            // I have no idea what broke doesn't work with the << operator. sleutenegger: not just Windows, but in general...???
            ARC_assert_stringstream << exceptionType <<  sfp.toString() << " " << message;
            throw(ARC_EXCEPTION_T(ARC_assert_stringstream.str()));
        }

        template<typename ARC_EXCEPTION_T>
        inline void ARC_throw_exception(std::string const & exceptionType, std::string const & function, std::string const & file,
                                        int line, std::string const & message)
        {
            ARC_throw_exception<ARC_EXCEPTION_T>(exceptionType, ARC::source_file_pos(function,file,line),message);
        }
    } // namespace ARC::detail

    template<typename ARC_EXCEPTION_T>
    inline void ARC_assert_throw(bool assert_condition, std::string message, ARC::source_file_pos sfp) {
        if(!assert_condition)
        {
            detail::ARC_throw_exception<ARC_EXCEPTION_T>("", sfp,message);
        }
    }
} // namespace ARC

#define ARC_THROW(exceptionType, message) {								\
    std::stringstream ARC_assert_stringstream;							\
    ARC_assert_stringstream << message;									\
    ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__, ARC_assert_stringstream.str()); \
  }


#define ARC_THROW_SFP(exceptionType, SourceFilePos, message){			\
    std::stringstream ARC_assert_stringstream;							\
    ARC_assert_stringstream << message;									\
    ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", SourceFilePos, ARC_assert_stringstream.str()); \
  }

#define ARC_ASSERT_TRUE(exceptionType, condition, message)				\
  if(!(condition))														\
    {																	\
      std::stringstream ARC_assert_stringstream;							         \
      ARC_assert_stringstream << "assert(" << #condition << ") failed: " << message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__, ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_FALSE(exceptionType, condition, message)				\
  if((condition))														\
    {																	\
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert( not " << #condition << ") failed: " << message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__, ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_GE_LT(exceptionType, value, lowerBound, upperBound, message) \
  if((value) < (lowerBound) || (value) >= (upperBound))							\
    {																        	\
      std::stringstream ARC_assert_stringstream;							    \
      ARC_assert_stringstream << "assert(" << #lowerBound << " <= " << #value << " < " << #upperBound << ") failed [" << (lowerBound) << " <= " << (value) << " < " << (upperBound) << "]: " << message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_LT(exceptionType, value, upperBound, message)			\
  if((value) >= (upperBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert(" << #value << " < " << #upperBound << ") failed [" << (value) << " < " << (upperBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_GE(exceptionType, value, lowerBound, message)			\
  if((value) < (lowerBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert(" << #value << " >= " << #lowerBound << ") failed [" << (value) << " >= " << (lowerBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_LE(exceptionType, value, upperBound, message)			\
  if((value) > (upperBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert(" << #value << " <= " << #upperBound << ") failed [" << (value) << " <= " << (upperBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_GT(exceptionType, value, lowerBound, message)			\
  if((value) <= (lowerBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert(" << #value << " > " << #lowerBound << ") failed [" << (value) << " > " << (lowerBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_EQ(exceptionType, value, testValue, message)			    \
  if((value) != (testValue))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert(" << #value << " == " << #testValue << ") failed [" << (value) << " == " << (testValue) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_NE(exceptionType, value, testValue, message)		     	\
  if((value) == (testValue))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "assert(" << #value << " != " << #testValue << ") failed [" << (value) << " != " << (testValue) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_NEAR(exceptionType, value, testValue, abs_error, message) \
  if(!(fabs((testValue) - (value)) <= fabs(abs_error)))						 \
    {																	     \
      std::stringstream ARC_assert_stringstream;							 \
      ARC_assert_stringstream << "assert(" << #value << " == " << #testValue << ") failed [" << (value) << " == " << (testValue) << " (" << fabs((testValue) - (value)) << " > " << fabs(abs_error) << ")]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#ifndef NDEBUG

#define ARC_THROW_DBG(exceptionType, message){							\
    std::stringstream ARC_assert_stringstream;							\
    ARC_assert_stringstream << message;									\
    ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__, ARC_assert_stringstream.str()); \
  }

#define ARC_ASSERT_TRUE_DBG(exceptionType, condition, message)			\
  if(!(condition))														\
    {																	\
      std::stringstream ARC_assert_stringstream;						\
      ARC_assert_stringstream << "debug assert(" << #condition << ") failed: " << message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__, ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_FALSE_DBG(exceptionType, condition, message)			\
  if((condition))														\
    {																	\
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert( not " << #condition << ") failed: " << message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__, ARC_assert_stringstream.str()); \
    }


#define ARC_ASSERT_DBG_RE( condition, message) ARC_ASSERT_DBG(std::runtime_error, condition, message)

#define ARC_ASSERT_GE_LT_DBG(exceptionType, value, lowerBound, upperBound, message) \
  if((value) < (lowerBound) || (value) >= (upperBound))							    \
    {																	            \
      std::stringstream ARC_assert_stringstream;							        \
      ARC_assert_stringstream << "debug assert(" << #lowerBound << " <= " << #value << " < " << #upperBound << ") failed [" << (lowerBound) << " <= " << (value) << " < " << (upperBound) << "]: " << message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_LT_DBG(exceptionType, value, upperBound, message)		\
  if((value) >= (upperBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert(" << #value << " < " << #upperBound << ") failed [" << (value) << " < " << (upperBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_GE_DBG(exceptionType, value, lowerBound, message)		\
  if((value) < (lowerBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert(" << #value << " >= " << #lowerBound << ") failed [" << (value) << " >= " << (lowerBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_LE_DBG(exceptionType, value, upperBound, message)		\
  if((value) > (upperBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert(" << #value << " <= " << #upperBound << ") failed [" << (value) << " <= " << (upperBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }

#define ARC_ASSERT_GT_DBG(exceptionType, value, lowerBound, message)		\
  if((value) <= (lowerBound))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert(" << #value << " > " << #lowerBound << ") failed [" << (value) << " > " << (lowerBound) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_EQ_DBG(exceptionType, value, testValue, message)		    \
  if((value) != (testValue))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert(" << #value << " == " << #testValue << ") failed [" << (value) << " == " << (testValue) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }


#define ARC_ASSERT_NE_DBG(exceptionType, value, testValue, message)		    \
  if((value) == (testValue))												\
    {																	    \
      std::stringstream ARC_assert_stringstream;							\
      ARC_assert_stringstream << "debug assert(" << #value << " != " << #testValue << ") failed [" << (value) << " != " << (testValue) << "]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }



#define ARC_ASSERT_NEAR_DBG(exceptionType, value, testValue, abs_error, message) \
  if(!(fabs((testValue) - (value)) <= fabs(abs_error)))						     \
    {																	         \
      std::stringstream ARC_assert_stringstream;						      	 \
      ARC_assert_stringstream << "debug assert(" << #value << " == " << #testValue << ") failed [" << (value) << " == " << (testValue) << " (" << fabs((testValue) - (value)) << " > " << fabs(abs_error) << ")]: " <<  message; \
      ARC::detail::ARC_throw_exception<exceptionType>("[" #exceptionType "] ", __FUNCTION__,__FILE__,__LINE__,ARC_assert_stringstream.str()); \
    }


#define ARC_OUT(X) std::cout << #X << ": " << (X) << std::endl

#else

#define ARC_OUT(X)
#define ARC_THROW_DBG(exceptionType, message)
#define ARC_ASSERT_TRUE_DBG(exceptionType, condition, message)
#define ARC_ASSERT_FALSE_DBG(exceptionType, condition, message)
#define ARC_ASSERT_GE_LT_DBG(exceptionType, value, lowerBound, upperBound, message)
#define ARC_ASSERT_LT_DBG(exceptionType, value, upperBound, message)
#define ARC_ASSERT_GT_DBG(exceptionType, value, lowerBound, message)
#define ARC_ASSERT_LE_DBG(exceptionType, value, upperBound, message)
#define ARC_ASSERT_GE_DBG(exceptionType, value, lowerBound, message)
#define ARC_ASSERT_NE_DBG(exceptionType, value, testValue, message)
#define ARC_ASSERT_EQ_DBG(exceptionType, value, testValue, message)
#define ARC_ASSERT_NEAR_DBG(exceptionType, value, testValue, abs_error, message)
#endif

#endif // ARC_ASSERT_MACROS_HPP

