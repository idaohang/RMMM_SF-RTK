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

#ifndef ARC_SOURCE_FILE_POS_HPP
#define ARC_SOURCE_FILE_POS_HPP

#include <string>
#include <iostream>
#include <sstream>
/// A class and macro that gives you the current file position.

/// \brief ARC Main namespace of this package.
namespace ARC {

class source_file_pos
{
public:
    std::string function;
    std::string file;
    int line;

    source_file_pos(std::string function, std::string file, int line) :
        function(function), file(file), line(line) {}

    operator std::string()
    {
        return toString();
    }

    std::string toString() const
    {
        std::stringstream s;
        s << file << ":" << line << ": " << function << "()";
        return s.str();
    }
  };

}// namespace ARC

inline std::ostream & operator<<(std::ostream & out, const ARC::source_file_pos & sfp)
{
    std::cout << sfp.file << ":" << sfp.line << ": " << sfp.function << "()";
    return out;
}
#define ARC_SOURCE_FILE_POS ARC::source_file_pos(__FUNCTION__,__FILE__,__LINE__)
#endif // ARC_SOURCE_FILE_POS_HPP
