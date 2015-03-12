#pragma once
#include <stdint.h>
#include <SDL.h>
#include <array>
#define GLEW_STATIC
#include <GL/glew.h>
#include <memory>

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

static std::array<uint8_t, 3> getNewColor(int value, int total = 10) {
	return HsvToRgb({ (uint8_t)(255.0f*(static_cast<float>(value) / total)), 128, 255 });
}
class SDLWrapper {
public:
	SDLWrapper(const char* name, int x, int y, size_t width, size_t height, size_t winScale = 1)
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

	~SDLWrapper()
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

#define printOpenGLError() printOglError(__FILE__, __LINE__)

int printOglError(char *file, int line)
{

	GLenum glErr;
	int    retCode = 0;

	glErr = glGetError();
	if (glErr != GL_NO_ERROR)
	{
		printf("glError in file %s @ line %d: %s\n",
			file, line, gluErrorString(glErr));
		retCode = 1;
	}
	return retCode;
}

class SDLWrapperGL {
public:
	SDLWrapperGL(const char* name, int x, int y, size_t width, size_t height, size_t winScale = 1, int drawMode = DRAW_MODE_IMAGE)
		: x(x), y(y), w(width), h(height), depth_(width*height), rgb_(3 * width*height), normals_(width*height), points_(3 * width*height), vbo_(5 * width*height), elements(width*height)
	{
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
		SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
		win = SDL_CreateWindow(name, x, y, (int)(width*winScale), (int)(height*winScale), SDL_WINDOW_OPENGL);
		context = SDL_GL_CreateContext(win);
		SDL_GL_MakeCurrent(win, context);
		switch (drawMode) {
		case DRAW_MODE_IMAGE:
			initGLState();
			break;
		case DRAW_MODE_POINTS:
			initGLPoints();
			break;
		case DRAW_MODE_QUADS:
			break;
		}
	}
	template <size_t subsamplePower>
	void setGrayScale(uint16_t *data) {
		static std::vector<uint8_t> tmp(4 * w*h);
		uint8_t *dest = (uint8_t *)tmp.data();
		for (int i = 0; i < w*h; i++)
			dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = (data[i] >> subsamplePower) > 255 ? 255 : (data[i] >> subsamplePower);

		setRGBA(dest);
	}

	void setGrayScale(uint8_t *data) {
		static std::vector<uint8_t> tmp(4 * w*h);
		uint8_t *dest = (uint8_t *)tmp.data();

		for (int i = 0; i < w*h; i++)
			dest[4 * i + 2] = dest[4 * i + 1] = dest[4 * i] = data[i];

		setRGBA(dest);
	}

	void setRGBA(uint8_t *data) {
		SDL_GL_MakeCurrent(win, context);
		glUseProgram(shaderProgram);

		GLuint tex;
		glGenTextures(1, &tex);
		glBindTexture(GL_TEXTURE_2D, tex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		// Draw a rectangle from the 2 triangles using 6 indices											
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		glDeleteTextures(1, &tex);

	}

	void setNormals(float *normals) {
		static std::vector<uint8_t> tmp(4 * w*h);
		uint8_t *dest = (uint8_t *)tmp.data();
		for (int i = 0; i < w*h; i++) {
			dest[4 * i] = static_cast<uint8_t>(128 + normals[3 * i] * 127);
			dest[4 * i + 1] = static_cast<uint8_t>(128 + normals[3 * i + 1] * 127);
			dest[4 * i + 2] = static_cast<uint8_t>(128 + normals[3 * i + 2] * 127);
		}
		setRGBA(dest);
	}
	void handleEvent(const SDL_Event & e){
		if (SDL_GetWindowFromID(e.window.windowID) == win) {
			//std::cout << "my win" << std::endl;
		}
	}
	void setPoints(float *points, float * normals, uint8_t * rgb) {																
		SDL_GL_MakeCurrent(win, context);																						
		glUseProgram(shaderProgram);																						  
		elements.clear();																						  				
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);																			  		
		glClear(GL_COLOR_BUFFER_BIT);
		int x = 0, y = 0;																									
		for (int i = 0; i < w*h; i++) {																						
			vbo_[5 * i] = points[3 * i];														
			vbo_[5 * i + 1] = points[3 * i+1];													
			vbo_[5 * i + 2] = points[3 * i+2];																				
			vbo_[5 * i + 3] = x++ / ((float)w);																				
			vbo_[5 * i + 4] = y / ((float)h);																		  			
			if (x == w) {																									  
				x = 0;																										  
				y++;																										  
			}																												  
			if (points[3 * i + 2] == 0)	{																						
				continue;																									  
			}																													
			elements.push_back(i);																							  
																															  
		}																													  
		glBindBuffer(GL_ARRAY_BUFFER, vbo);																					  
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*vbo_.size(), vbo_.data(), GL_STREAM_DRAW);									
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);																			  
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*elements.size(), elements.data(), GL_STREAM_DRAW);
		GLuint tex;	
		float vFov = 1.0, ar = 1.333, nearclip = 500, farclip = 650000;
		const auto yf = 1 / std::tan(vFov / 2);
		const auto xf = yf / ar;
		const auto dz = nearclip - farclip;
		float cameraMat[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1/80.0, 0 };
		glUniformMatrix4fv(cameraAttrib, 1, false, cameraMat);

		{																													  
			static std::vector<uint8_t> tmp(3 * w*h);																		  
			uint8_t *dest = (uint8_t *)tmp.data();																			  
			for (int i = 0; i < w*h; i++) {																					  
				dest[3 * i] = static_cast<uint8_t>(128 + normals[3 * i] * 127);												  
				dest[3 * i + 1] = static_cast<uint8_t>(128 + normals[3 * i + 1] * 127);										  
				dest[3 * i + 2] = static_cast<uint8_t>(128 + normals[3 * i + 2] * 127);										  
			}																												  
			glGenTextures(1, &tex);																							  
			glBindTexture(GL_TEXTURE_2D, tex);																				  
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);												  
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, tmp.data());							  
  
		}																													  
		glDrawElements(GL_POINTS, elements.size(), GL_UNSIGNED_INT, 0);
		glDeleteTextures(1, &tex);																						  		

		printOpenGLError();
	}																														  
	void Draw() {
		SDL_GL_MakeCurrent(win, context);
		SDL_GL_SwapWindow(win);
	}
	~SDLWrapperGL()
	{
		deinitGLState();
		SDL_DestroyWindow(win);
	}
	enum {
		DRAW_MODE_POINTS,
		DRAW_MODE_QUADS,
		DRAW_MODE_IMAGE
	};
private:
	std::vector<uint16_t> depth_;
	std::vector<float> points_;
	std::vector<float> normals_;
	std::vector<uint8_t> rgb_;

	std::vector<float> vbo_;
	std::vector<GLuint> elements;				
	size_t w, h;
	int x, y;
	SDL_Window * win;
	SDL_GLContext context;
	GLuint vao, vbo, ebo;
	GLuint vertexShader, fragmentShader, shaderProgram;
	GLuint posAttrib, texAttrib,transAttrib,cameraAttrib;
	void deinitGLState()
	{
		SDL_GL_MakeCurrent(win, context);
		glDeleteProgram(shaderProgram);
		glDeleteShader(fragmentShader);
		glDeleteShader(vertexShader);

		glDeleteBuffers(1, &ebo);
		glDeleteBuffers(1, &vbo);

		glDeleteVertexArrays(1, &vao);
	}
	void initGLState()
	{
		glewExperimental = GL_TRUE;
		GLenum glewErr = glewInit();
		if (GLEW_OK != glewErr)
		{
			std::cerr << "Failed to initialize GLEW." << std::endl;
			std::cerr << glewGetErrorString(glewErr) << std::endl;
		}

		const GLchar* vertexSource =
			"#version 430 core\n"
			"in vec2 position;"
			"in vec2 texcoord;"
			"out vec2 Texcoord;"
			"void main() {"
			"   Texcoord = texcoord;"
			"   gl_Position = vec4(position, 0.0, 1.0);"
			"}";
		const GLchar* fragmentSource =
			"#version 430 core\n"
			"in vec2 Texcoord;"
			"out vec4 outColor;"
			"uniform sampler2D tex;"
			"void main() {"
			"   outColor = texture(tex, Texcoord);"
			"}";

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		// Create a Vertex Buffer Object and copy the vertex data to it
		glGenBuffers(1, &vbo);

		static GLfloat vertices[] = {
			//  Position   Color             Texcoords
			-1.f, 1.f, 0.0f, 0.0f, // Top-left
			1.f, 1.f, 1.0f, 0.0f, // Top-right
			1.f, -1.f, 1.0f, 1.0f, // Bottom-right
			-1.f, -1.f, 0.0f, 1.0f  // Bottom-left
		};

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

		// Create an element array																								 
		glGenBuffers(1, &ebo);

		static GLuint elements[] = {
			0, 1, 2,
			2, 3, 0
		};

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

		// Create and compile the vertex shader																					 
		vertexShader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertexShader, 1, &vertexSource, NULL);
		glCompileShader(vertexShader);

		// Create and compile the fragment shader																				 
		fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
		glCompileShader(fragmentShader);

		// Link the vertex and fragment shader into a shader program															 
		shaderProgram = glCreateProgram();
		glAttachShader(shaderProgram, vertexShader);
		glAttachShader(shaderProgram, fragmentShader);
		glBindFragDataLocation(shaderProgram, 0, "outColor");
		glLinkProgram(shaderProgram);

		// Specify the layout of the vertex data																				 
		posAttrib = glGetAttribLocation(shaderProgram, "position");
		glEnableVertexAttribArray(posAttrib);
		glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

		texAttrib = glGetAttribLocation(shaderProgram, "texcoord");
		glEnableVertexAttribArray(texAttrib);
		glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));
	}
	void initGLPoints()
	{
		glewExperimental = GL_TRUE;
		GLenum glewErr = glewInit();
		if (GLEW_OK != glewErr)
		{
			std::cerr << "Failed to initialize GLEW." << std::endl;
			std::cerr << glewGetErrorString(glewErr) << std::endl;
		}

		const GLchar* vertexSource =
			"#version 430 core\n"
			"in vec3 position;"
			"in vec2 texcoord;"
			"out vec2 Texcoord;"
			"uniform mat4 trans;"
			"uniform mat4 camera;"
			"void main() {"
			"   Texcoord = texcoord;"
			"   float scaledZ = position.z/250;"
			"   gl_PointSize = 7.0;"//clamp(scaledZ,0.25,8);"
			"   vec2 xy = 2.0*position.xy/position.z;"
			"   xy.x = xy.x*3.0/4.0;"
			"   xy.y = -xy.y;"
			"   gl_Position = camera*vec4(position,1.0);" //vec4(xy, 1.0, 1.0)
			"}";
		const GLchar* fragmentSource =
			"#version 430 core\n"
			"in vec2 Texcoord;"
			"out vec4 outColor;"
			"uniform sampler2D tex;"
			"void main() {"
			//"   outColor = vec4(1.0,0.0,1.0,1.0);"
			"   outColor = texture(tex, Texcoord);"
			"}";
		glEnable(GL_PROGRAM_POINT_SIZE);
		//glPointSize(3);
		glGenVertexArrays(1, &vao);																								 
		glBindVertexArray(vao);																									 
																																 
		// Create a Vertex Buffer Object and copy the vertex data to it															 
		glGenBuffers(1, &vbo);																									 
		// Create an element array																								 
		glGenBuffers(1, &ebo);																									 
		glBindBuffer(GL_ARRAY_BUFFER, vbo);																					  	
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*vbo_.size(), vbo_.data(), GL_STREAM_DRAW);									
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);																			  	
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*elements.size(), elements.data(), GL_STREAM_DRAW);
		vertexShader = glCreateShader(GL_VERTEX_SHADER);																		 
		glShaderSource(vertexShader, 1, &vertexSource, NULL);																	 
		glCompileShader(vertexShader);																							 
		GLint success = 0;
		glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
		if (success == GL_FALSE) {
			GLint maxLength = 0;
			glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH, &maxLength);

			// The maxLength includes the NULL character
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(vertexShader, maxLength, &maxLength, &errorLog[0]);
			std::cout << errorLog.data() << std::endl;
		}
																																 
		// Create and compile the fragment shader																				 
		fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);																	 
		glShaderSource(fragmentShader, 1, &fragmentSource, NULL);																 
		glCompileShader(fragmentShader);																						 
		success = 0;
		glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
		if (success == GL_FALSE) {
			GLint maxLength = 0;
			glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &maxLength);

			// The maxLength includes the NULL character
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(fragmentShader, maxLength, &maxLength, &errorLog[0]);
			std::cout << errorLog.data() << std::endl;
		}
																																 
		// Link the vertex and fragment shader into a shader program															 
		shaderProgram = glCreateProgram();																						 
		glAttachShader(shaderProgram, vertexShader);																			 
		glAttachShader(shaderProgram, fragmentShader);																			 
		glBindFragDataLocation(shaderProgram, 0, "outColor");																	 
		glLinkProgram(shaderProgram);																							 
		glUseProgram(shaderProgram);																							 

		// Specify the layout of the vertex data																				 
		posAttrib = glGetAttribLocation(shaderProgram, "position");																 
		glEnableVertexAttribArray(posAttrib);																					 
		glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0);										 
		printOpenGLError();																										
		texAttrib = glGetAttribLocation(shaderProgram, "texcoord");																 
		glEnableVertexAttribArray(texAttrib);																					 
		glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));			
		transAttrib = glGetUniformLocation(shaderProgram, "trans");
		cameraAttrib = glGetUniformLocation(shaderProgram, "camera");
	}
};


