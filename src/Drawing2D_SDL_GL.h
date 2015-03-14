#pragma once
#include <stdint.h>
#include <SDL.h>
#include <array>
#define GLEW_STATIC
#include <GL/glew.h>
#include <memory>
#include <vector>

class Draw2DImage {
public:
	Draw2DImage(const char* name, int x, int y, size_t width, size_t height, size_t winScale = 1)
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
		initGLState();
		
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
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (GLsizei)w, (GLsizei)h, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
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
	void Draw() {
		SDL_GL_MakeCurrent(win, context);
		SDL_GL_SwapWindow(win);
	}
	~Draw2DImage()
	{
		deinitGLState();
		SDL_DestroyWindow(win);
	}
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
	GLuint posAttrib, texAttrib, transAttrib, cameraAttrib;
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
};


