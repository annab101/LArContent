/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the charge splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ChargeSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ChargeSplittingAlgorithm::ChargeSplittingAlgorithm() :
    m_chargeThreshold(1.0f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChargeSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, CartesianVector &splitPosition) const
{
    bool foundSplit(false);
    bool foundClusterSplit(false);
    CartesianVector splitClusterPosition(0.f, 0.f, 0.f);
    float distToCluster = 9999.;

    std::cout << "Looking for charge split!" << std::endl;

    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    std::cout << orderedCaloHitList.size() << std::endl;

    float prevHitEnergy = -1;

    for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
    {   
        if(foundClusterSplit == true){
            continue;
        }

        for (const CaloHit *const pCaloHit : *layerEntry.second)
        {
            if(foundClusterSplit == true){
                continue;
            }
           float hitEnergy = pCaloHit->GetInputEnergy();

           if(prevHitEnergy > 0 && std::abs(hitEnergy - prevHitEnergy) > 1000){
                foundClusterSplit = true;
                std::cout << "split cluster!!" << std::endl;
                std::cout << "prev energy: " << prevHitEnergy << std::endl;
                std::cout << "current energy: " << hitEnergy << std::endl;
                splitClusterPosition = pCaloHit->GetPositionVector();
                std::cout << pCaloHit->GetPositionVector().GetX() << " "  << pCaloHit->GetPositionVector().GetY() << " "  << pCaloHit->GetPositionVector().GetZ() << " "  << std::endl;
           }

           //std::cout << prevHitEnergy << std::endl;
           //std::cout << "difference: " << std::abs(hitEnergy - prevHitEnergy) << std::endl;
           prevHitEnergy = hitEnergy;
           //std::cout << pCaloHit->GetPositionVector().GetX() << " "  << pCaloHit->GetPositionVector().GetY() << " "  << pCaloHit->GetPositionVector().GetZ() << " "  << std::endl;
           
        }
    }

    if(foundClusterSplit){

        std::cout << "finding layer split" << std::endl;

        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());
        for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
        {
            const int iLayer(iter->first);

            const float rL(slidingFitResult.GetL(iLayer));
            
            CartesianVector centralPosition(0.f, 0.f, 0.f);


            if (STATUS_CODE_SUCCESS != slidingFitResult.GetGlobalFitPosition(rL, centralPosition))
            {
                continue;
            }

            const CartesianVector clusterToLayer(centralPosition - splitClusterPosition);
            const float distanceSquared(clusterToLayer.GetMagnitudeSquared());
            //std::cout << "dist squared: " << distanceSquared << std::endl;
            //std::cout << centralPosition.GetX() << " " << centralPosition.GetY() << " " << centralPosition.GetZ() << std::endl;

                if (distanceSquared < distToCluster)
                {
                    distToCluster = distanceSquared;
                    splitPosition = centralPosition;
                    foundSplit = true;
                }

        }

        std::cout << "Split at " << splitPosition.GetX() << " " << splitPosition.GetY() << " " << splitPosition.GetZ() << std::endl;

    }



    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChargeSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChargeThreshold", m_chargeThreshold));

    return TwoDSlidingFitSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
