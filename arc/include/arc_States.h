//
// Created by sujinglan on 7/7/17.
//

#ifndef ARC_ARC_STATES_H
#define ARC_ARC_STATES_H
#include "arc.h"
/**
 * @class ARC_States
 * @brief the SRTK states for particel filter
 *
 * This state has the following parameters:
 * @li <b>xb</b>    the X-position of the base station (ECEF/m)
 * @li <b>yb</b>    the Y-Position of the base station (ECEF/m)
 * @li <b>zb</b>    the Z-Position of the base station (ECEF/m)
 * @li <b>xu</b>    the X-position of the rover station (ECEF/m)
 * @li <b>yu</b>    the Y-position of the rover station (ECEF/m)
 * @li <b>zu</b>    the Z-position of the rover station (ECEF/m)
 * @li <b>dtb</b>   the clock drift of the base station (m)
 * @li <b>dtu</b>   the clock drift of the rover station (m)
 * @li <b>Ni</b>    the ambiguity of GPS and BDS ,every satelite have a ambiguity (cycle)
 * @li <b>Tb</b>    the trosphere delay of the base station
 * @li <b>Tu</b>    the trosphere delay of the rover station
 * @author sujinglan
 */
class ARC_States{
public:
    /**
     * Constructor
     */
    ARC_States();
    /**
     * Destructor
     */
    virtual ~ARC_States();

    /**
     * Opereator
     */
    ARC_States operator*(double factor) const;
    ARC_States& operator+=(const ARC_States& other);

    /// \brief set the position of base station in ECEF
    void setBasePosX(double BasePosX);
    void setBasePosY(double BasePosY);
    void setBasePosZ(double BasePosZ);

    /// \brief set the position of rover station in ECEF
    void setRoverPosX(double RoverX);
    void setRoverPosY(double RoverY);
    void setRoverPosZ(double RoverZ);

    /// \brief get the postion of base station in ECEF
    double getBasePosX() const;
    double getBasePosY() const;
    double getBasePosZ() const;

    /// \brief get the postion of rover station in ECEF
    double getRoverPosX() const;
    double getRoverPosY() const;
    double getRoverPosZ() const;

    /// \brief set the clock drift value of base and rover station (m)
    void setBaseClk(double BaseClk);
    void setRoverClk(double RoverClk);

    /// \brief get the clock drift value of base and rover station (m)
    double getRoverClk() const;
    double getBaseClk() const ;

    /// \brief set the trosphere delay of base and rover station
    void setBaseTrop(double BaseTrop);
    void setRoverTrop(double RoverTrop);

    /// \brief get the trosphere delay of base and rover station
    double getBaseTrop() const ;
    double getRoverTrop() const ;

    /// \brief set the ambguity of every satelite for single-difference on base-rover station
    void setAmb(double *N);

    /// \brief get the ambguity
    double* getAmb() const;

private:
    /// \brief the states for base station position
    double m_BaseX;
    double m_BaseY;
    double m_BaseZ;
    /// \brief the states for rover station postion
    double m_RoverX;
    double m_RoverY;
    double m_RoverZ;
    /// \brief the clock drift of base and rover station
    double m_BaseClk;
    double m_RoverClk;
    /// \brief the trosphere delay of base station
    double m_BaseTrop;
    double m_RoverTrop;
    /// \brief the ambiguity value of every satelite
    double m_N[MAXSAT];
};
#endif //ARC_ARC_STATES_H
