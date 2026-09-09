/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeRatioSplittingAlgorithm.h
 *
 *  @brief  Header file for the overshoot splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHARGERATIO_SPLITTING_ALGORITHM_H
#define LAR_CHARGERATIO_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitMultiSplitAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ChargeRatioSplittingAlgorithm class
 */
class ChargeRatioSplittingAlgorithm : public TwoDSlidingFitMultiSplitAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ChargeRatioSplittingAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FindBestSplitPositions(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterSplittingMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minClusterLength;      ///<
    float m_maxIntersectDisplacement;
    float m_minSplitDisplacement;
};

} // namespace lar_content

#endif // #ifndef LAR_OVERSHOOT_SPLITTING_ALGORITHM_H
