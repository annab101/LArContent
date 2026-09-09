/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeRatioSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the charge ratio splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeRatioSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ChargeRatioSplittingAlgorithm::ChargeRatioSplittingAlgorithm() :
    TwoDSlidingFitMultiSplitAlgorithm(),
    m_minClusterLength(5.f),
    m_maxIntersectDisplacement(1.5f),
    m_minSplitDisplacement(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChargeRatioSplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
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

void ChargeRatioSplittingAlgorithm::FindBestSplitPositions(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterSplittingMap) const
{
    std::cout << "running charge split based on ratio" << std::endl;
    ClusterList clusterList;
    for (const auto &mapEntry : slidingFitResultMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : clusterList)
    {
        bool foundClusterSplit(false);
        CartesianVector splitClusterPosition(0.f, 0.f, 0.f);
        CartesianVector prevSplitClusterPosition(-9999.f, -9999.f, -9999.f);
        float prevSplitSize = 0.;
        CaloHitList caloHitList;
        clusterList.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        int caloHitList_size = caloHitList.size();

        if(caloHitList_size < 7)
        {
            continue;
        }
        
        for(int i=3; i<caloHitList_size-3; i++){
            auto it_prev3 = std::next(caloHitList.begin(), i-3);
            auto it_prev2 = std::next(caloHitList.begin(), i-2);
            auto it_prev1 = std::next(caloHitList.begin(), i-1);
            auto it_next3 = std::next(caloHitList.begin(), i+3);
            auto it_next2 = std::next(caloHitList.begin(), i+2);
            auto it_next1 = std::next(caloHitList.begin(), i+1);
            auto it = std::next(caloHitList.begin(), i);
            float prev3HitAv = (float)((*it_prev3)->GetInputEnergy() + (*it_prev2)->GetInputEnergy() + (*it_prev1)->GetInputEnergy())/3;
            float next3HitAv = (float)((*it_next3)->GetInputEnergy() + (*it_next2)->GetInputEnergy() + (*it_next1)->GetInputEnergy())/3;

            std::cout << prev3HitAv << " | " << next3HitAv << std::endl;

            if((prev3HitAv > 1.5*next3HitAv) || (1.5*prev3HitAv < next3HitAv)){
                foundClusterSplit = true;
                std::cout << "=====================================" << std::endl;
                std::cout << "prev average: " << prev3HitAv << std::endl;
                std::cout << "current average: " << next3HitAv << std::endl;
                splitClusterPosition = (*it)->GetPositionVector();
                std::cout << (*it)->GetPositionVector().GetX() << " "  << (*it)->GetPositionVector().GetY() << " "  << (*it)->GetPositionVector().GetZ() << " "  << std::endl;
            }

            const CartesianVector splitDistanceVect(splitClusterPosition - prevSplitClusterPosition);
            const float distanceSquared(splitDistanceVect.GetMagnitudeSquared());
            float prevSplitDistance = sqrt(distanceSquared);
            float splitSize = prev3HitAv/next3HitAv;

            if(foundClusterSplit && prevSplitDistance > 1.0)
            {
                std::cout << "split cluster!!" << std::endl;
                clusterSplittingMap[pCluster].push_back(splitClusterPosition);
            }
            else if(foundClusterSplit && prevSplitDistance < 1.0){
                //choose best split for split position within 1cm
                if(std::abs(splitSize - 1) > std::abs(prevSplitSize - 1)){
                    std::cout << "split cluster!!" << std::endl;
                    clusterSplittingMap[pCluster].push_back(splitClusterPosition);
                }
                else{
                    std::cout << "split cluster!!" << std::endl;
                    clusterSplittingMap[pCluster].push_back(prevSplitClusterPosition);
                }
            }

            foundClusterSplit = false;
            prevSplitClusterPosition = splitClusterPosition;
            prevSplitSize = splitSize;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChargeRatioSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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
