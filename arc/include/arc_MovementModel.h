//
// Created by sujinglan on 7/7/17.
//

#ifndef ARC_ARC_MOVEMENTMODEL_H
#define ARC_ARC_MOVEMENTMODEL_H
#include <libPF/MovementModel.h>
#include "arc_States.h"

/**
 * @class ARC_MovementModel
 *
 * @brief single frequency rtk position movement model propagates class
 *
 * This movement model propagates a car state according to the translation and
 * rotation speed.
 *
 * @author sujinglan
 */
class ARC_MovementModel : public libPF::MovementModel<ARC_States>
{
    /**
     * Constructor
     */
    ARC_MovementModel();

    /**
     * Destructor
     */
    ~ARC_MovementModel();

    /**
     * The drift method propagates the car using its speed.
     * @param state Pointer to the state that has to be manipulated.
     */
    void drift(ARC_States& state, double dt) const;

    /**
     * The diffusion consists of a very small gaussian jitter on the
     * state's variable.
     * @param state Pointer to the state that has to be manipulated.
     */
    void diffuse(ARC_States& state, double dt) const;

    /**
     * @param d new standard deviation for the diffusion of x,y,z
     *          of base and rover station position
     */
    void setBaesXStdDev(double d) {m_Std_BaseX=d;}
    void setBaseYstdDev(double d) {m_Std_BaseY=d;}
    void setBaseZStdDev(double d) {m_Std_BaseZ=d;}
    void setRoverXStdDev(double d) {m_Std_RoverX=d;}
    void setRoverYStdDev(double d) {m_Std_RoverY=d;}
    void setRoverZStdDev(double d) {m_Std_RoverZ=d;}
    
    /**
     * @param d new standard deviation for the diffusion of clock drift
     *          of base and rover station
     */
    void setBaseClkStdDev(double d) {m_Std_Base_Clk=d;}
    void setRoverClkStdDev(double d) {m_Std_Rover_Clk=d;}
    /**
     * @param d new standard deviation for the diffusion of trosphere delay
     *          of base and rover station
     */
    void setBaseTrpStdDev(double d) {m_Std_Base_Trp=d;}
    void setRoverTrpStdDev(double d) {m_Std_Rover_Trp=d;}
    /**
     * @param d new standard deviation for the diffusion of ambguity
     *          of single difference ambguity on base and rover station
     */
    void setNStdDev(double *d) {
        for (size_t i=0;i<MAXSAT;i++) m_Std_N[i]=d[i];
    }

private:
    /// \brief standard deviation for the base and rover position
    double m_Std_BaseX,m_Std_BaseY,m_Std_BaseZ;
    double m_Std_RoverX,m_Std_RoverY,m_Std_RoverZ;
    /// \brief standard deviation for the base and rover clock drift
    double m_Std_Base_Clk,m_Std_Rover_Clk;
    /// \brief standard deviation for the base and rover trosphere delay
    double m_Std_Base_Trp,m_Std_Rover_Trp;
    /// \brief standard deviation for ambguity
    double m_Std_N[MAXSAT];

    /// \brief Stores the random number generator
    libPF::RandomNumberGenerationStrategy* m_RNG;

    /**
     * Private Method
     */
    double SQR(double x) { return x*x;}

};
#endif //ARC_ARC_MOVEMENTMODEL_H
