/**
 *  @file   larpandoracontent/LArMonitoring/MCHitMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MCHitMonitoringAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArObjects/LArEventTopology.h"


#include <sstream>

using namespace pandora;

namespace lar_content
{

MCHitMonitoringAlgorithm::MCHitMonitoringAlgorithm() :
    m_caloHitListName(""),
    m_mcListName("Input"),
    m_visualise(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCHitMonitoringAlgorithm::~MCHitMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCHitMonitoringAlgorithm::Run()
{

    this->BuildMCHitMap();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCHitMonitoringAlgorithm::BuildMCHitMap()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcListName, pMCParticleList));

    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_SUCCESS;

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    m_mcHitsMap.clear();
    for (const MCParticle *const pMC : *pMCParticleList)
        m_mcHitsMap[pMC] = CaloHitList();

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            m_mcHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }
    }

    MCHitsMap hitAboveThresholdMap, otherHitsMap;
    for (const auto &[pMC, hits] : m_mcHitsMap)
    {
        for (const CaloHit *const pCaloHit : hits)
        {
            bool above_threshold{false};
            if (pCaloHit->GetInputEnergy() > 1000.)
            {
                above_threshold = true;
                hitAboveThresholdMap[pMC].emplace_back(pCaloHit);
            }        
            if (!above_threshold)
            {
                otherHitsMap[pMC].emplace_back(pCaloHit);
            }
        }
    }

    for (const auto &[pMC, hits] : hitAboveThresholdMap)
    {
        if (hits.empty())
            continue;

        if (m_visualise)
        {
            std::ostringstream oss;
            oss << pMC->GetParticleId() << " (" << hits.size() << ")";
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, oss.str(), AUTOITER));
        }

        MCParticleWeightMap mcParticleWeightMap;

        for (const CaloHit *const pCaloHit : hits)
        {
            const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

            MCParticleVector mcParticleVector;
            for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap) mcParticleVector.push_back(mapEntry.first);
            std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

            std::cout << "Found hit belonging to " << pMC->GetParticleId() << " with input energy "
                    << pCaloHit->GetInputEnergy() << std::endl;
            for (const MCParticle *const pMCParticle : mcParticleVector)
            {
                const float weight(hitMCParticleWeightMap.at(pMCParticle));
                mcParticleWeightMap[pMCParticle] += weight;
                std::cout << "MC particle: " << pMCParticle->GetParticleId() << ", weight: " << weight << std::endl;
            }

        }
    }
    for (const auto &[pMC, hits] : otherHitsMap)
    {
        if (hits.empty())
            continue;

        if (m_visualise)
        {
            std::ostringstream oss;
            oss << pMC->GetParticleId() << " (" << hits.size() << ")";
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, oss.str(), BLACK));
        }
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    /*for (const auto &[pMC, hits] : m_mcHitsMap)
    {
        const MCParticleList &parentList{pMC->GetParentList()};
        std::string desc{""};
        if (parentList.size() == 1)
        {
            std::ostringstream oss;
            oss << parentList.front();
            desc += std::to_string(parentList.front()->GetParticleId()) + " (" + oss.str() + ")" + " -> ";
            if (parentList.front()->GetParticleId() == 111)
            {
                std::cout << "MCParticle: PDG " << pMC->GetParticleId() << " with parent PDG 111 has " << hits.size() << " hits." << std::endl;
            }
        }
        desc += std::to_string(pMC->GetParticleId());

        if (m_visualise && !hits.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, desc, AUTOITER));
        }
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }*/

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCHitMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCListName", m_mcListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content