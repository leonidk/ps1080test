#pragma once
#include <stdint.h>
#include <SDL.h>
#include <array>

class SDLWrapper {
public:
	SDLWrapper(const char* name, int x, int y, size_t width, size_t height,size_t winScale = 1)
		: x(x), y(y), w(width), h(height)
	{
		win = SDL_CreateWindow(name, x, y, (int)(width*winScale), (int)(height*winScale), 0);
		ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
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
		for (int i = 0; i <  sur->h * sur->w; i++) {
			dest[4 * i] = data[3*i];
			dest[4 * i+1] = data[3*i+1];
			dest[4 * i+2] = data[3*i+2];

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
	static std::array<uint8_t, 3> getNewColor(int value, int total = 10) {
		return HsvToRgb({ (uint8_t)(255.0f*(static_cast<float>(value) / total)), 128, 255 });
	}
	~SDLWrapper()
	{
		SDL_DestroyRenderer(ren);
		SDL_DestroyWindow(win);
		SDL_FreeSurface(sur);
	}
private:
	static std::array<uint8_t, 3> HsvToRgb(std::array<uint8_t, 3> hsv)
	{
		std::array<uint8_t, 3> rgb;
		unsigned char region, remainder, p, q, t;

		if (hsv[1] == 0)
		{
			rgb[0] = hsv[2];
			rgb[1] = hsv[2];
			rgb[2] = hsv[2];
			return rgb;
		}

		region = hsv[0] / 43;
		remainder = (hsv[0] - (region * 43)) * 6;

		p = (hsv[2] * (255 - hsv[1])) >> 8;
		q = (hsv[2] * (255 - ((hsv[1] * remainder) >> 8))) >> 8;
		t = (hsv[2] * (255 - ((hsv[1] * (255 - remainder)) >> 8))) >> 8;

		switch (region)
		{
		case 0:
			rgb[0] = hsv[2]; rgb[1] = t; rgb[2] = p;
			break;
		case 1:
			rgb[0] = q; rgb[1] = hsv[2]; rgb[2] = p;
			break;
		case 2:
			rgb[0] = p; rgb[1] = hsv[2]; rgb[2] = t;
			break;
		case 3:
			rgb[0] = p; rgb[1] = q; rgb[2] = hsv[2];
			break;
		case 4:
			rgb[0] = t; rgb[1] = p; rgb[2] = hsv[2];
			break;
		default:
			rgb[0] = hsv[2]; rgb[1] = p; rgb[2] = q;
			break;
		}

		return rgb;
	}

	size_t w, h;
	int x, y;
	SDL_Window * win;
	SDL_Renderer * ren;
	SDL_Surface * sur;
};


