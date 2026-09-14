#pragma once

#include <string>

namespace ui {

// Loads languages/<name>.json (UTF-8, flat "key": "value" object).
// Missing keys keep the built-in English fallbacks.
bool loadLanguage(char const *name);

wchar_t const *tr(char const *key);
std::wstring trf(char const *key);
std::wstring trf(char const *key, std::wstring const &a0);
std::wstring trf(char const *key, std::wstring const &a0, std::wstring const &a1);

} // namespace ui
