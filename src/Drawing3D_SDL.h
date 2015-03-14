#pragma once
#include <stdint.h>
#include <SDL.h>
#include <array>
#define GLEW_STATIC
#include <GL/glew.h>
#include <memory>
#include "linalg.h"


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

class Draw3DImage {
public:
	Draw3DImage(const char* name, int x, int y, size_t width, size_t height, size_t winScale = 1, int drawMode = DRAW_MODE_POINTS)
		: x(x), y(y), w(width), h(height), depth_(width*height), rgb_(3 * width*height), normals_(width*height), points_(3 * width*height), vbo_(5 * width*height), elements(width*height)
	{
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
		SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
		SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
		win = SDL_CreateWindow(name, x, y, (int)(width*winScale), (int)(height*winScale), SDL_WINDOW_OPENGL);
		context = SDL_GL_CreateContext(win);
		SDL_GL_MakeCurrent(win, context);
		switch (drawMode) {
		case DRAW_MODE_POINTS:
			initGLPoints();
			break;
		case DRAW_MODE_QUADS:
			break;
		}
	}
	void handleEvent(const SDL_Event & e){
		auto MOVE_SCALE = 35.0f;
		if (SDL_GetWindowFromID(e.window.windowID) == win) {
			if (e.type == SDL_MOUSEMOTION) {
				if (e.motion.state & SDL_BUTTON_RMASK) {
					position[0] += MOVE_SCALE*e.motion.xrel;
					position[1] -= MOVE_SCALE*e.motion.yrel;
				} else if (e.motion.state & SDL_BUTTON_LMASK) {
					using namespace linalg;
					using namespace linalg::aliases;

					float4 xr(std::sin(e.motion.yrel / 180.0), 0, 0, std::cos(e.motion.yrel / 180.0));
					float4 yr(0, std::sin(e.motion.xrel / 180.0), 0, std::cos(e.motion.xrel / 180.0));
					auto combined = qmul(qmul(orientation,xr),yr);
					orientation = combined;
				} else if (e.motion.state & SDL_BUTTON_MMASK) {
					position[2] -= MOVE_SCALE*e.motion.yrel;
				}
			} else if (e.type == SDL_MOUSEWHEEL) {
				position[2] += MOVE_SCALE * e.wheel.y;
			}
		}
		if (e.type == SDL_KEYDOWN) {
			switch (e.key.keysym.scancode)
			{
			case SDL_SCANCODE_W:
				position += MOVE_SCALE* linalg::qzdir(orientation);
				break;
			case SDL_SCANCODE_A:
				position += MOVE_SCALE*linalg::qxdir(orientation);
				break;
			case SDL_SCANCODE_S:
				position -= MOVE_SCALE* linalg::qzdir(orientation);
				break;
			case SDL_SCANCODE_D:
				position -= MOVE_SCALE*linalg::qxdir(orientation);
				break;
			default:
				break;
			}
		
		}
	}
	void setPoints(float *points, float * normals, uint8_t * rgb, bool quads=true) {																
		auto dThresh = 0.05;
		SDL_GL_MakeCurrent(win, context);																						
		glUseProgram(shaderProgram);																						  
		elements.clear();																						  				
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);																			  		
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_CULL_FACE);
		if (quads) {
			auto cnt = 0;
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					auto i = y*w + x;
					vbo_[5 * i] = points[3 * i];
					vbo_[5 * i + 1] = points[3 * i + 1];
					vbo_[5 * i + 2] = points[3 * i + 2];
					vbo_[5 * i + 3] = x/((float)w);
					vbo_[5 * i + 4] = y/((float)h);

				}
			}
			for (int y = 0; y < h-1; y++) {
				for (int x = 0; x < w - 1; x++) {
					auto i = y*w + x;
					auto tL = points[3 * i + 2];
					auto bR = points[3 * ((y + 1)*w + x + 1) + 2];

					if (tL < 100 || bR < 100)
						continue;
					auto tR = points[3 * (i + 1) + 2];
					auto mean = (tL + tR + bR)*(1.0f / 3.0f);
					if (tR > 100 && abs(mean - tL) < dThresh*mean) {
						elements.push_back(y*w + x);
						elements.push_back((y + 1)*w + x + 1);
						elements.push_back(y*w + x + 1);
					}
					auto bL = points[3 * ((y + 1)*w + x) + 2];
					auto mean2 = (bR + bL + tL)*(1.0f / 3.0f);
					if (bL > 100 && abs(mean2 - tL) < dThresh*mean2) {
						elements.push_back((y + 1)*w + x + 1);
						elements.push_back(y*w + x);
						elements.push_back((y + 1)*w + x);
					}
				}
			}
		}
		else {
			int x = 0, y = 0;
			for (int i = 0; i < w*h; i++) {
				vbo_[5 * i] = points[3 * i];
				vbo_[5 * i + 1] = points[3 * i + 1];
				vbo_[5 * i + 2] = points[3 * i + 2];
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
		}
		glBindBuffer(GL_ARRAY_BUFFER, vbo);																					  
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*vbo_.size(), vbo_.data(), GL_STREAM_DRAW);									
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);																			  
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*elements.size(), elements.data(), GL_STREAM_DRAW);
		GLuint tex;	
		float vFov = 0.82, ar = 1.333, nearclip = 100, farclip = 100000;
		const auto yf = 1 / std::tan(vFov / 2);
		const auto xf = yf / ar;
		const auto dz = farclip - nearclip;
		
		float cameraMat[16] = { xf, 0, 0, 0, 0, yf, 0, 0, 0, 0, -(farclip + nearclip) / dz, -2*nearclip*farclip/dz, 0, 0, -1, 0 };
		glUniformMatrix4fv(cameraAttrib, 1, GL_TRUE, cameraMat);
		linalg::aliases::float4x4 a;
		auto qx = qxdir(orientation);
		auto qy = qydir(orientation);
		auto qz = qzdir(orientation);

		a[0] = { qx[0], qx[1], qx[2], 0 };
		a[1] = { qy[0], qy[1], qy[2], 0 };
		a[2] = { qz[0], qz[1], qz[2], 0 };
		a[3] = { position[0], position[1], position[2], 1 };		
		glUniformMatrix4fv(transAttrib, 1, GL_FALSE, (float*)&a);
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
		glEnable(GL_DEPTH_TEST);
		if (quads) {
			glDrawElements(GL_TRIANGLES, elements.size(), GL_UNSIGNED_INT, 0);
		}
		else {
			glDrawElements(GL_POINTS, elements.size(), GL_UNSIGNED_INT, 0);
		}
		glDisable(GL_DEPTH_TEST);
		glDeleteTextures(1, &tex);																						  		

		printOpenGLError();
	}																														  
	void Draw() {
		SDL_GL_MakeCurrent(win, context);
		SDL_GL_SwapWindow(win);
	}
	~Draw3DImage()
	{
		deinitGLState();
		SDL_DestroyWindow(win);
	}
	enum {
		DRAW_MODE_POINTS,
		DRAW_MODE_QUADS
	};
private:
	std::vector<uint16_t> depth_;
	std::vector<float> points_;
	std::vector<float> normals_;
	std::vector<uint8_t> rgb_;

	std::vector<float> vbo_;
	std::vector<GLuint> elements;		
	linalg::aliases::float3 position = { 0, 0, 0 };
	linalg::aliases::float4 orientation = { 0, 0, 0, 1 };

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
	void initGLPoints()
	{
		glewExperimental = GL_TRUE;
		GLenum glewErr = glewInit();

		const GLchar* vertexSource =
			"#version 430 core\n"
			"in vec3 position;"
			"in vec2 texcoord;"
			"out vec2 Texcoord;"
			"uniform mat4 trans;"
			"uniform mat4 camera;"
			"void main() {"
			"   Texcoord = texcoord;"
			//"   Texcoord = position.yz;"
			"  vec4 cameraCoord =trans*vec4(position.x,-position.yz,1.0); "
			"   gl_PointSize = clamp(-cameraCoord.z/400,4,4);"
			" vec4 outp = camera*cameraCoord;"
			"   gl_Position = outp;"//vec4(outp.xy,1.0,1.0);" //vec4(xy, 1.0, 1.0)
			"}";
		const GLchar* fragmentSource =
			"#version 430 core\n"
			"in vec2 Texcoord;"
			"out vec4 outColor;"
			"uniform sampler2D tex;"
			"void main() {"
			//"   outColor = vec4(1000.0/Texcoord.y,1000.0/Texcoord.y,1000.0/Texcoord.y,1.0);"
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


