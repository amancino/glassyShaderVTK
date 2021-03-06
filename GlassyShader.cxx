#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkNamedColors.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkShaderProgram.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkTriangleMeshPointNormals.h>
#include <vtkVersion.h>


#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtksys/SystemTools.hxx>
#include <vtkUniforms.h>
#include <vtkOpenGLCamera.h>
#include <vtkOpenGLTexture.h>
#include <vtkTextureObject.h>
#include <vtkJPEGReader.h>
#include <vtkImageFlip.h>

#if VTK_VERSION_NUMBER >= 89000000000ULL
#define USE_SHADER_PROPERTIES 1
#include <vtkShaderProperty.h>
#endif

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName);


template <typename T>
void ToFloat(T* in, float* out, int noOfComponents)
{
  for (int i = 0; i < noOfComponents; ++i)
  {
    out[i] = static_cast<float>(in[i]);
  }
}

/*
MC: Model coordinates: local
WC: World coordinates
VC: View coordinates: from the eyes (camera)
DC: Device coordinates: display

Normal matrix:
The normal matrix is the matrix which preserves vertex normals under an affine transform.
If you do the math, this turns out to be the inverse transpose of the modelview matrix.
Note that you can make the code much more efficient whenever the transform, i.e. the modelview,
happens to be orthogonal; this would imply that the inverse transpose is just the identity
transform, which then makes the normal matrix the same as the modelview matrix.

The following link describes coordinates systems in openGL:
https://learnopengl.com/Getting-started/Coordinate-Systems

Model matrix:
Maps from MC (model/local coordinates) to WC.
I believe this is internal to VTK, since actors exist in WC.

View matrix:
Maps from WC to VC. Commonly refered to as "ModelView matrix".

Projection matrix:
Maps from WC to DC.

*/



void SetCameraShaderParameters(vtkShaderProgram* prog, vtkRenderer* ren, vtkOpenGLCamera* cam)
{
  vtkMatrix4x4* glTransformMatrix;
  vtkMatrix4x4* modelViewMatrix;
  vtkMatrix3x3* normalMatrix;
  vtkMatrix4x4* projectionMatrix;
  // modelViewMatrix = WCVCMatrix: WC to VC (world to view) transformation
  // normalMatrix: transforms normals (vectors) from WC to VC
  // projectionMatrix = VCDCMatrix: VC to DC (view to display) transformation
  // glTransformMatrix = WCDCMatrix: WC to DC (world to display) transformation
  cam->GetKeyMatrices(ren, modelViewMatrix, normalMatrix, projectionMatrix, glTransformMatrix);

  vtkNew<vtkMatrix4x4> InverseProjectionMat;
  vtkNew<vtkMatrix4x4> InverseModelViewMat;
  vtkNew<vtkMatrix4x4> InverseVolumeMat;

  InverseProjectionMat->DeepCopy(projectionMatrix);
  InverseProjectionMat->Invert();
 // prog->SetUniformMatrix("in_projectionMatrix", projectionMatrix);
 // prog->SetUniformMatrix("in_inverseProjectionMatrix", InverseProjectionMat.GetPointer());

  InverseModelViewMat->DeepCopy(modelViewMatrix);
  InverseModelViewMat->Invert();
 // prog->SetUniformMatrix("in_modelViewMatrix", modelViewMatrix);
 // prog->SetUniformMatrix("in_inverseModelViewMatrix", InverseModelViewMat.GetPointer());

  float fvalue3[3];
  if (cam->GetParallelProjection())
  {
    double dir[4];
    cam->GetDirectionOfProjection(dir);
    ToFloat(dir, fvalue3, 3);
    //prog->SetUniform3fv("in_projectionDirection", 1, &fvalue3);
  }

  // camera position in WC
  ToFloat(cam->GetPosition(), fvalue3, 3);
  prog->SetUniform3fv("u_camera", 1, &fvalue3);
  //cout << "u_camera: " << fvalue3[0] << ", " << fvalue3[1] << ", " << fvalue3[2] << endl;
}


class vtkShaderCallback : public vtkCommand
{
public:
  static vtkShaderCallback* New()
  {
    return new vtkShaderCallback;
  }
  vtkRenderer* Renderer;
  vtkOpenGLTexture* InputTexture;
  void Execute(vtkObject*, unsigned long, void* calldata) override
  {
    vtkShaderProgram* program = reinterpret_cast<vtkShaderProgram*>(calldata);
    if (program)
    {
      // update camera position
      vtkOpenGLCamera* cam = vtkOpenGLCamera::SafeDownCast(Renderer->GetActiveCamera());
      SetCameraShaderParameters(program,Renderer,cam);

      InputTexture->GetTextureObject()->Activate();
      int unit = InputTexture->GetTextureUnit();
      program->SetUniformi("u_cubemap", InputTexture->GetTextureUnit());

    }
  }

  vtkShaderCallback()
  {
    this->Renderer = nullptr;
    this->InputTexture = nullptr;
  }
};



//----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  auto polyData = ReadPolyData("/home/axel/Desktop/BSCI/patients/bronchi.stl");
  auto colors = vtkSmartPointer<vtkNamedColors>::New();

  auto actor = vtkSmartPointer<vtkActor>::New();
  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  auto mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
  renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

  auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(640, 480);
  renderWindow->AddRenderer(renderer);
  renderer->AddActor(actor);

  auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);

  auto triangles = vtkSmartPointer<vtkTriangleFilter>::New();
  triangles->SetInputData(polyData);
  triangles->Update();

  auto normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetInputData(triangles->GetOutput());
  normals->SetAutoOrientNormals(1);
  normals->Update();

  auto norms = vtkSmartPointer<vtkTriangleMeshPointNormals>::New();
  norms->SetInputData(triangles->GetOutput());
  norms->Update();

  mapper->SetInputData(normals->GetOutput());
  mapper->ScalarVisibilityOff();

  actor->SetMapper(mapper);

  vtkShaderProperty* sp = actor->GetShaderProperty();




  actor->GetProperty()->SetAmbientColor(0.2, 0.2, 0.2);
  actor->GetProperty()->SetDiffuseColor(1.0, 1.0, 1.0);
  actor->GetProperty()->SetSpecularColor(1.0, 1.0, 1.0);
  actor->GetProperty()->SetSpecular(0.5);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetAmbient(0.1);
  //actor->GetProperty()->SetSpecularPower(100.0);
  actor->GetProperty()->SetOpacity(0.5);

  // set a new opacity uniform
  float opacity = actor->GetProperty()->GetOpacity();

  // shader works properly only with with vtk opacity in 1.0
  actor->GetProperty()->SetOpacity(1.0);


  // read vertex shader (modified version, adapted to vtk)
  std::ifstream vertexShader("/home/axel/work/tutorial/VTK/GlassyShader/glass.vtkVert.glsl");
  std::ostringstream vertexCode;
  vertexCode << vertexShader.rdbuf();

  // read fragment shader (modified version, adapted to vtk)
  std::ifstream fragmentShader("/home/axel/work/tutorial/VTK/GlassyShader/glass.vtkFrag.glsl");
  std::ostringstream fragmentCode;
  fragmentCode << fragmentShader.rdbuf();

  std::cout << "vertex code: " << std::endl
            << vertexCode.str() << std::endl;

  std::cout << "fragment code: " << std::endl
            << fragmentCode.str() << std::endl;

  // Modify the vertex shader to pass the position of the vertex
  sp->AddVertexShaderReplacement(
      "//VTK::Normal::Dec",  // replace the normal block
      true,                  // before the standard replacements
      "//VTK::Normal::Dec\n" // we still want the default
      //"uniform mat4 u_viewProjectionMatrix;\n"  // MCDCMatrix? or VCDCMatrix?
      //"uniform mat4 u_modelMatrix;\n"           // MCVCMatrix? should be MCWCMatrix
      //"uniform mat3 u_normalMatrix;\n"  // normalMatrix
      "uniform vec4 u_camera;\n"          // camera position in WC

      //"in vec4 a_vertex;\n"             // vertexMC (vertex position in MC)
      //"in vec3 a_normal;\n"             // normalMC (vertex normal in MC)
      "out vec3 v_reflection;\n"
      "out vec3 v_refraction;\n"
      "out float v_fresnel;\n"
      // Indices of refraction
      "const float Air = 1.0;\n"
      "const float Glass = 1.51714;\n"
      // Air to glass ratio of the indices of refraction (Eta)
      "const float Eta = Air / Glass;\n"
      // see http://en.wikipedia.org/wiki/Refractive_index Reflectivity
      "const float R0 = ((Air - Glass) * (Air - Glass)) / ((Air + Glass) * (Air + Glass));\n",
      false // only do it once
  );

  sp->AddVertexShaderReplacement(
      "//VTK::Normal::Impl",  // replace the normal block
      true,                   // before the standard replacements
      "//VTK::Normal::Impl\n" // we still want the default
      "vec4 vertex = MCVCMatrix*vertexMC;\n"  // vertex in VC (should it be WC?)

      //"vec3 incident = normalize(vec3(vertex-u_camera));\n" // bug: camera is in WC and vertex in VC
      "vec4 cameraVC = MCVCMatrix*u_camera;\n"
      "vec3 incident = normalize(vec3(vertex-cameraVC));\n" // in VC

      // Assume incoming normal is normalized.
      //"vec3 normal = normalMatrix*normalMC;\n"    // normal in VC
      // normal in VC is already computed in "normalVCVSOutput"

      "v_refraction = refract(incident, normalVCVSOutput, Eta);\n"
      "v_reflection = reflect(incident, normalVCVSOutput);\n"

      // see http://en.wikipedia.org/wiki/Schlick%27s_approximation
      "v_fresnel = R0 + (1.0 - R0) * pow((1.0 - dot(-incident, normalVCVSOutput)), 5.0);\n",

     // "gl_Position = MCDCMatrix*vertexMC;\n",   // vertex position in DC
      false // only do it once
  );

  // Define varying and uniforms for the fragment shader here
  sp->AddFragmentShaderReplacement(
      "//VTK::Normal::Dec",  // replace the normal block
      true,                  // before the standard replacements
      "//VTK::Normal::Dec\n" // we still want the default
      "uniform samplerCube u_cubemap;\n"
      //"uniform float opacity;\n"  // declared using "SetUniformf"
      "in vec3 v_refraction;\n"
      "in vec3 v_reflection;\n"
      "in float v_fresnel;\n"
      "out vec4 fragColor;\n",
      false // only do it once
  );

  sp->AddFragmentShaderReplacement(
        "//VTK::Color::Dec",
        true,
        "\n",
        //"//VTK::Color::Dec\n",  // ambient, opacity, specular, etc.
        false
  );

  sp->AddFragmentShaderReplacement(
        "//VTK::Color::Impl",
        true,
        "\n",
        //"//VTK::Color::Impl\n",
        false
  );

  sp->AddFragmentShaderReplacement(
      "//VTK::Light::Impl",  // replace the light block
      true,                 // after the standard replacements
      "  vec4 refractionColor = texture(u_cubemap, normalize(v_refraction));\n"
      "  vec4 reflectionColor = texture(u_cubemap, normalize(v_reflection));\n"
      "\n"
      "  fragColor = mix(refractionColor, reflectionColor, v_fresnel);\n"
      "  fragOutput0 = fragColor;\n"
      "  fragOutput0.a = opacity;\n",
      //"  ambientColor = fragColor.rgb;\n"
      //"  opacity = fragColor.a;\n"
      //"//VTK::Light::Impl\nn", // we still want the default calc
      //"  fragOutput0= fragColor;\n"
      //"  fragOutput0.a = opacityUniform;\n",
      true // only do it once
  );

  auto myCallback = vtkSmartPointer<vtkShaderCallback>::New();
  myCallback->Renderer = renderer;
  mapper->AddObserver(vtkCommand::UpdateShaderEvent, myCallback);

  // add cube texture map
  vtkNew<vtkOpenGLTexture> texture;
  texture->CubeMapOn();
  texture->UseSRGBColorSpaceOn();

  /*
  std::string pathSkybox[6] = { "../Data/skybox/posx.jpg", "../Data/skybox/negx.jpg",
    "../Data/skybox/posy.jpg", "../Data/skybox/negy.jpg", "../Data/skybox/posz.jpg",
    "../Data/skybox/negz.jpg" };
*/
  std::string pathSkybox[6] = { "../Data/skybox2/front.jpg", "../Data/skybox2/right.jpg",
    "../Data/skybox2/top.jpg", "../Data/skybox2/bottom.jpg", "../Data/skybox2/back.jpg",
    "../Data/skybox2/left.jpg" };

  for (int i = 0; i < 6; i++)
  {
    vtkNew<vtkJPEGReader> jpg;

    jpg->SetFileName(pathSkybox[i].c_str());
    vtkNew<vtkImageFlip> flip;
    flip->SetInputConnection(jpg->GetOutputPort());
    flip->SetFilteredAxis(1); // flip y axis
    texture->SetInputConnection(i, flip->GetOutputPort());
  }
  myCallback->InputTexture = texture;
  myCallback->Renderer->SetEnvironmentTexture(myCallback->InputTexture);
  myCallback->Renderer->UseImageBasedLightingOn();


  sp->GetFragmentCustomUniforms()->SetUniformf("opacity",opacity);

  renderWindow->Render();
  renderer->GetActiveCamera()->SetPosition(-.3, 0, .08);
  renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
  renderer->GetActiveCamera()->SetViewUp(.26, 0.0, .96);
  renderer->ResetCamera();
  renderer->GetActiveCamera()->Zoom(0.5);
  renderWindow->Render();
  renderWindow->SetWindowName("GlassyShader");
  renderWindow->Render();
  interactor->Start();
  return EXIT_SUCCESS;
}

vtkSmartPointer<vtkPolyData> ReadPolyData(const char* fileName)
{
  vtkSmartPointer<vtkPolyData> polyData;
  std::string extension =
      vtksys::SystemTools::GetFilenameExtension(std::string(fileName));
  if (extension == ".ply")
  {
    auto reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".vtp")
  {
    auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".obj")
  {
    auto reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".stl")
  {
    auto reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".vtk")
  {
    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else if (extension == ".g")
  {
    auto reader = vtkSmartPointer<vtkBYUReader>::New();
    reader->SetGeometryFileName(fileName);
    reader->Update();
    polyData = reader->GetOutput();
  }
  else
  {
    auto source = vtkSmartPointer<vtkSphereSource>::New();
    source->SetPhiResolution(25);
    source->SetThetaResolution(25);
    source->Update();
    polyData = source->GetOutput();
  }
  return polyData;
}
