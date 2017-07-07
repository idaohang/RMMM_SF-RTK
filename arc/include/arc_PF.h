//
// Created by sujinglan on 7/7/17.
//

#ifndef ARC_ARC_PF_H
#define ARC_ARC_PF_H
/**
 * @brief ARC-SRTK Particle Filter Some Const-Variable and Marco
 *
 * @author sujinglan
 */
#define ARC_PF_BASEPOS_STD                 (1.0)           /// initial standard deviation of base station position (m)
#define ARC_PF_ROVERPOS_STD                (100.0)         /// initial standard deviation of rover station position (m)
#define ARC_PF_BASECLK_STD                 (200.0)         /// initial standard deviation of base station clock drift (m)
#define ARC_PF_ROVERCLK_STD                (200.0)         /// initial standard deviation of rover station clock drift (m)
#define ARC_PF_BASETROP_STD                (5.0)           /// initial standard deviation of base station trosphere delay (m)
#define ARC_PF_ROVERTROP_STD               (5.0)           /// initial standard deviation of rover station trosphere delay (m)
#define ARC_PF_AMB_STD                     (10.0)          /// initial standard deviation of single-difference ambguity (cycle)

#endif //ARC_ARC_PF_H
