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
        //float prevHitEnergy = -1;
        bool foundClusterSplit(false);
        CartesianVector splitClusterPosition(0.f, 0.f, 0.f);
        CartesianVector prevSplitClusterPosition(-9999.f, -9999.f, -9999.f);
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
            float prev3HitSum = (*it_prev3)->GetInputEnergy() + (*it_prev2)->GetInputEnergy() + (*it_prev1)->GetInputEnergy();
            float next3HitSum = (*it_next3)->GetInputEnergy() + (*it_next2)->GetInputEnergy() + (*it_next1)->GetInputEnergy();

            if((prev3HitSum > 4000 && next3HitSum < 3000) || (prev3HitSum < 3000 && next3HitSum > 4000)){
                foundClusterSplit = true;
                std::cout << "=====================================" << std::endl;
                std::cout << "prev average: " << prev3HitSum << std::endl;
                std::cout << "current average: " << next3HitSum << std::endl;
                splitClusterPosition = (*it)->GetPositionVector();
                std::cout << (*it)->GetPositionVector().GetX() << " "  << (*it)->GetPositionVector().GetY() << " "  << (*it)->GetPositionVector().GetZ() << " "  << std::endl;
            }

            const CartesianVector splitDistanceVect(splitClusterPosition - prevSplitClusterPosition);
            const float distanceSquared(splitDistanceVect.GetMagnitudeSquared());
            float prevSplitDistance = sqrt(distanceSquared);

            if(foundClusterSplit && prevSplitDistance > .9)
            {
                std::cout << "split cluster!!" << std::endl;
                clusterSplittingMap[pCluster].push_back(splitClusterPosition);
            }
            
            foundClusterSplit = false;
            prevSplitClusterPosition = splitClusterPosition;
        }

        /*float prevHitEnergy = -1.;

        for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
        {   
            for (const CaloHit *const pCaloHit : *layerEntry.second)
            {
                bool foundClusterSplit(false);
                CartesianVector splitClusterPosition(0.f, 0.f, 0.f);
                float hitEnergy = pCaloHit->GetInputEnergy();
                float hitDif = std::abs(hitEnergy - prevHitEnergy);

                if(prevHitEnergy > 0 && hitDif > 500){
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
        }*/
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
