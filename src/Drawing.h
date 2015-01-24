#pragma once
#include <stdint.h>
#include <SDL.h>

void drawGrayScale(SDL_Surface* sur, uint16_t * data) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
		dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = (data[i] >> 1) > 255 ? 255 : (data[i] >> 1);
}

void drawGrayScale(SDL_Surface* sur, uint8_t * data) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
		dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = data[i];
}

void drawNormals(SDL_Surface* sur, float* normals) {
	uint8_t* dest = (uint8_t*)sur->pixels;
	for (int i = 0; i < sur->h*sur->w; i++)
	{
		dest[4 * i] = 128 + normals[3 * i] * 127;
		dest[4 * i + 1] = 128 + normals[3 * i + 1] * 127;
		dest[4 * i + 2] = 128 + normals[3 * i + 2] * 127;

	}
}


