#pragma once
#include <stdint.h>
#include <SDL.h>
#include <array>

class Draw2DImage {
public:
	Draw2DImage(const char* name, int x, int y, size_t width, size_t height, size_t winScale = 1)
		: x(x), y(y), w(width), h(height)
	{
		win = SDL_CreateWindow(name, x, y, (int)(width*winScale), (int)(height*winScale), 0);
		ren = SDL_CreateRenderer(win, -1, 0);
		sur = SDL_CreateRGBSurface(0, (int)width, (int)height, 32, 0, 0, 0, 0);

	}
	template <size_t subsamplePower>
	void setGrayScale(uint16_t *data) {
		uint8_t *dest = (uint8_t *)sur->pixels;
		for (int i = 0; i < sur->h * sur->w; i++)
			dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = (data[i] >> subsamplePower) > 255 ? 255 : (data[i] >> subsamplePower);
	}

	void setGrayScale(uint8_t *data) {
		uint8_t *dest = (uint8_t *)sur->pixels;
		for (int i = 0; i < sur->h * sur->w; i++)
			dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = data[i];
	}

	void setRGB(uint8_t *data) {
		uint8_t *dest = (uint8_t *)sur->pixels;
		for (int i = 0; i < sur->h * sur->w; i++) {
			dest[4 * i] = data[3 * i];
			dest[4 * i + 1] = data[3 * i + 1];
			dest[4 * i + 2] = data[3 * i + 2];

		}
	}

	void setNormals(float *normals) {
		uint8_t *dest = (uint8_t *)sur->pixels;
		for (int i = 0; i < sur->h * sur->w; i++) {
			dest[4 * i] = static_cast<uint8_t>(128 + normals[3 * i] * 127);
			dest[4 * i + 1] = static_cast<uint8_t>(128 + normals[3 * i + 1] * 127);
			dest[4 * i + 2] = static_cast<uint8_t>(128 + normals[3 * i + 2] * 127);
		}
	}
	void Draw()
	{
		SDL_Texture *tex = SDL_CreateTextureFromSurface(ren, sur);
		SDL_RenderClear(ren);
		SDL_RenderCopy(ren, tex, NULL, NULL);
		SDL_RenderPresent(ren);
		SDL_DestroyTexture(tex);
	}

	~Draw2DImage()
	{
		SDL_DestroyRenderer(ren);
		SDL_DestroyWindow(win);
		SDL_FreeSurface(sur);
	}
private:
	size_t w, h;
	int x, y;
	SDL_Window * win;
	SDL_Renderer * ren;
	SDL_Surface * sur;
};