/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeMultiSplittingAlgorithm.h
 *
 *  @brief  Header file for the overshoot splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHARGEMULTI_SPLITTING_ALGORITHM_H
#define LAR_CHARGEMULTI_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitMultiSplitAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ChargeMultiSplittingAlgorithm class
 */
class ChargeMultiSplittingAlgorithm : public TwoDSlidingFitMultiSplitAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ChargeMultiSplittingAlgorithm();

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
