#pragma once

#include "chemistry/SaceSimulationCompiler.h"
#include "substance/RuntimeSubstance.h"

#include <cstddef>
#include <cstdint>

// SACE Phase 18: session-local generated runtime material cache.
// Authoritative generated identity remains SaceCatalog / SaceGeneratedRecord.
// This registry maps a compact RuntimeSubstanceRef (never-reused handle)
// to SaceRecordId + compiled simulation profile.
//
// Lifetime is tied to SaceCatalog: SaceCatalog::clear() also resets this
// registry. Handles are never reused in-process, so stale refs cannot resolve
// to a later record that reused a SaceRecordId.
//
// Registration is not spawnability. Generated records stay spawnable = false.
// Do not call from physics ticks. Built-ins are not stored here.

struct GeneratedRuntimeMaterial {
    RuntimeSubstanceRef runtimeRef{};
    SaceRecordId sourceRecord = kSaceRecordNone;
    SaceCompiledSimulationProfile profile{};
    uint32_t displayOrdinal = 0;
    bool valid = false;
};

void resetGeneratedRuntimeMaterials();
std::size_t generatedRuntimeMaterialCount();

bool registerGeneratedRuntimeMaterial(SaceRecordId recordId, RuntimeSubstanceRef &out);

// Lookup already-registered generated catalog id. Does not register.
RuntimeSubstanceRef runtimeGenerated(SaceRecordId recordId);

SaceRecordId runtimeGeneratedId(RuntimeSubstanceRef ref);

RuntimeSubstanceRef runtimeRefFromSaceSubstanceRef(SaceSubstanceRef const &ref);

GeneratedRuntimeMaterial const *generatedRuntimeMaterial(RuntimeSubstanceRef ref);

SaceCompiledLiquidProfile const *runtimeLiquidProfile(RuntimeSubstanceRef ref);
SaceCompiledGasProfile const *runtimeGasProfile(RuntimeSubstanceRef ref);
SaceCompiledPhaseProfile const *runtimeLiquidGasPhaseProfile(RuntimeSubstanceRef ref);

// Writes a display label. Not canonical identity. Not a lookup key.
bool runtimeFormatDisplayName(RuntimeSubstanceRef ref, char *out, int cap);

// Headless: --sace-runtime-registry-diag -> misc/sace_runtime_registry_diag.tsv
void runSaceRuntimeRegistryDiagnostics();
