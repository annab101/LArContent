/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeMultiSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the overshoot splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeMultiSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ChargeMultiSplittingAlgorithm::ChargeMultiSplittingAlgorithm() :
    TwoDSlidingFitMultiSplitAlgorithm(),
    m_minClusterLength(5.f),
    m_maxIntersectDisplacement(1.5f),
    m_minSplitDisplacement(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChargeMultiSplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChargeMultiSplittingAlgorithm::FindBestSplitPositions(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterSplittingMap) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : slidingFitResultMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : clusterList)
    {
        float prevHitEnergy = -1;

        for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
        {   
            for (const CaloHit *const pCaloHit : *layerEntry.second)
            {
                bool foundClusterSplit(false);
                CartesianVector splitClusterPosition(0.f, 0.f, 0.f);
                float hitEnergy = pCaloHit->GetInputEnergy();
                float hitDif = std::abs(hitEnergy - prevHitEnergy);

                if(prevHitEnergy > 0 && hitDif > 200){
                        foundClusterSplit = true;
                        std::cout << "split cluster!!" << std::endl;
                        std::cout << "prev energy: " << prevHitEnergy << std::endl;
                        std::cout << "current energy: " << hitEnergy << std::endl;
                        splitClusterPosition = pCaloHit->GetPositionVector();
                        std::cout << pCaloHit->GetPositionVector().GetX() << " "  << pCaloHit->GetPositionVector().GetY() << " "  << pCaloHit->GetPositionVector().GetZ() << " "  << std::endl;
                }

                prevHitEnergy = hitEnergy;
                    
                if(foundClusterSplit)
                {
                    clusterSplittingMap[pCluster].push_back(splitClusterPosition);
                }
                
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChargeMultiSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxIntersectDisplacement", m_maxIntersectDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSplitDisplacement", m_minSplitDisplacement));

    return TwoDSlidingFitMultiSplitAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
