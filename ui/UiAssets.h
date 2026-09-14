#pragma once

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

namespace ui {

void initAssets();
void shutdownAssets();

// Draws a category glyph into the left side of `slot`. Returns the width
// reserved for the icon (0 if that PNG is missing).
int drawCategoryIcon(HDC dc, RECT const &slot, int categoryIndex, bool dimmed);

} // namespace ui
