/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeSplittingAlgorithm.h
 *
 *  @brief  Header file for the kink splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHARGE_SPLITTING_ALGORITHM_H
#define LAR_CHARGE_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ChargeSplittingAlgorithm class
 */
class ChargeSplittingAlgorithm : public TwoDSlidingFitSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ChargeSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Use sliding linear fit to identify the best split position
     *
     *  @param  slidingFitResult the input sliding fit result
     *  @param  splitPosition the best split position
     *
     *  @return pandora::StatusCode
     */
    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, pandora::CartesianVector &splitPosition) const;

    float m_chargeThreshold;      ///<=
};

} // namespace lar_content

#endif // #ifndef LAR_KINK_SPLITTING_ALGORITHM_H
