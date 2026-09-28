#pragma once

#include "chemistry/SaceCatalog.h"

#include <cstdint>
#include <type_traits>

// SACE Phase 18: compact engine-facing substance handle.
// Built-in SubstanceId and generated SaceRecordId stay separate types.
// This packed 32-bit value is for future hot liquid/gas component storage.
// It is not canonical chemistry identity and is not a SubstanceId.
//
// Encoding (raw == 0 is none):
//   bits 31-30: kind  01 = built-in, 10 = generated
//   bits 29-0 : payload
//     built-in payload = SubstanceId (1 .. SUBSTANCE_COUNT-1)
//     generated payload = process-lifetime runtime handle (never a SaceRecordId)
// Generated SaceRecordId is stored only in the runtime registry entry.
// Handles are never reused after reset, so a stale ref cannot alias a new record.

enum class RuntimeSubstanceKind : uint8_t {
    None = 0,
    BuiltIn,
    Generated
};

struct RuntimeSubstanceRef {
    uint32_t raw = 0;
};

constexpr uint32_t kRuntimeKindShift = 30;
constexpr uint32_t kRuntimeKindBuiltIn = 1u;
constexpr uint32_t kRuntimeKindGenerated = 2u;
constexpr uint32_t kRuntimePayloadMask = 0x3FFFFFFFu;
constexpr uint32_t kRuntimeGeneratedHandleMax = kRuntimePayloadMask;

static_assert(std::is_trivially_copyable_v<RuntimeSubstanceRef>,
    "RuntimeSubstanceRef must be trivially copyable for component arrays");
static_assert(sizeof(RuntimeSubstanceRef) == 4,
    "RuntimeSubstanceRef is a fixed 32-bit packed handle");

inline bool operator==(RuntimeSubstanceRef a, RuntimeSubstanceRef b) {
    return a.raw == b.raw;
}
inline bool operator!=(RuntimeSubstanceRef a, RuntimeSubstanceRef b) {
    return a.raw != b.raw;
}

inline RuntimeSubstanceKind runtimeSubstanceKind(RuntimeSubstanceRef ref) {
    if (ref.raw == 0)
        return RuntimeSubstanceKind::None;
    uint32_t k = ref.raw >> kRuntimeKindShift;
    if (k == kRuntimeKindBuiltIn)
        return RuntimeSubstanceKind::BuiltIn;
    if (k == kRuntimeKindGenerated)
        return RuntimeSubstanceKind::Generated;
    return RuntimeSubstanceKind::None;
}

inline uint32_t runtimeSubstancePayload(RuntimeSubstanceRef ref) {
    return ref.raw & kRuntimePayloadMask;
}

inline RuntimeSubstanceRef runtimeNone() {
    return {};
}

inline RuntimeSubstanceRef runtimeBuiltIn(SubstanceId id) {
    if (id == SUBSTANCE_NONE || id >= SUBSTANCE_COUNT)
        return {};
    RuntimeSubstanceRef ref;
    ref.raw = (kRuntimeKindBuiltIn << kRuntimeKindShift) | static_cast<uint32_t>(id);
    return ref;
}

// Pack a never-reused generated runtime handle. This is not SaceRecordId.
inline RuntimeSubstanceRef runtimeGeneratedHandle(uint32_t handle) {
    if (handle == 0 || handle > kRuntimeGeneratedHandleMax)
        return {};
    RuntimeSubstanceRef ref;
    ref.raw = (kRuntimeKindGenerated << kRuntimeKindShift) | handle;
    return ref;
}

inline bool runtimeSubstanceIsNone(RuntimeSubstanceRef ref) {
    return runtimeSubstanceKind(ref) == RuntimeSubstanceKind::None;
}

inline bool runtimeSubstanceIsBuiltIn(RuntimeSubstanceRef ref) {
    if (runtimeSubstanceKind(ref) != RuntimeSubstanceKind::BuiltIn)
        return false;
    uint32_t v = runtimeSubstancePayload(ref);
    return v > 0 && v < static_cast<uint32_t>(SUBSTANCE_COUNT);
}

// Representation-valid generated encoding. May still be stale after registry reset.
inline bool runtimeSubstanceIsGenerated(RuntimeSubstanceRef ref) {
    return runtimeSubstanceKind(ref) == RuntimeSubstanceKind::Generated
        && runtimeSubstancePayload(ref) != 0;
}

inline bool runtimeSubstanceValid(RuntimeSubstanceRef ref) {
    return runtimeSubstanceIsBuiltIn(ref) || runtimeSubstanceIsGenerated(ref);
}

inline SubstanceId runtimeBuiltinId(RuntimeSubstanceRef ref) {
    if (!runtimeSubstanceIsBuiltIn(ref))
        return SUBSTANCE_NONE;
    return static_cast<SubstanceId>(runtimeSubstancePayload(ref));
}

inline uint32_t runtimeGeneratedHandleValue(RuntimeSubstanceRef ref) {
    if (!runtimeSubstanceIsGenerated(ref))
        return 0;
    return runtimeSubstancePayload(ref);
}
