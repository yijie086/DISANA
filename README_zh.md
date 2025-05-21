# **DISANA**

---

## **目录结构**

```
DISANA/
├── macros/
│   ├── Run101.C
│   ├── main101.C
├── source/
│   ├── Cuts/
│   │   ├── EventCut.h
│   │   ├── EventCut.cxx
│   ├── ParticleInformation/
│   │   ├── RECParticle.h
│   │   ├── RECParticle.cxx
│   ├── Math/
│   │   ├── RECParticleKinematic.h
│   │   ├── RECParticleKinematic.cxx
│   ├── core/
│   │   ├── FilesInPath.h
│   │   ├── FilesInPath.cxx
│   ├── DrawHist/
│       ├── DrawAndSave.h
│       ├── DrawAndSave.cxx
```

---

## **`macros`说明**

#### **1. macros/Run101.C**

- **功能**：
  - 这是主分析脚本，负责加载 HIPO 数据文件、应用筛选条件、计算粒子运动学变量、生成直方图并保存结果。
- **主要逻辑**：
  1. 调用 `GetHipoFilesInPath` 获取指定目录下的所有 `.hipo` 文件 (可遍历嵌套文件夹中所有`.hipo`文件)。
  2. 使用 `ROOT::RDataFrame` 加载数据。
  3. 定义筛选条件（如电子、质子和光子的筛选条件可自定义，细节参见`EventCut.cxx`和`EventCut.h`）。
  4. 计算粒子的运动学变量（如 `theta`, `phi`, `p`, 并且添加到 `ROOT::RDF::RNode` 上）。
  5. 调用`DrawAndSave.cxx`内的函数生成并保存直方图为`.png`和`.root`。
  6. 将筛选后的数据snapshot到 ROOT 文件中。

---

#### **2. macros/main101.C**

- **功能**：
  - 这是程序的入口文件，负责解析命令行参数并调用 `Run101` 函数。
- **主要逻辑**：
  1. 检查命令行参数是否提供了 HIPO 数据文件路径。
  2. 调用 `Run101` 函数执行数据分析。

---

## `source` 说明 

#### **1. Cuts/EventCut.h 和 `Cuts/EventCut.cxx`**

- **功能**：
  - 定义了一个通用的事件筛选器 EventCut
  - 用于根据粒子的电荷、动量、顶点位置等条件筛选事件。

- **主要逻辑**：
  1. 提供设置筛选条件的函数（如 `SetChargeCut`、`SetMomentumCut` 等）。
  2. 支持与 `ROOT::RDataFrame` 的无缝集成。

- **依赖**：
  - 无直接依赖，但与数据结构配合使用。

---

#### **2. `ParticleInformation/RECParticle.h` 和 `ParticleInformation/RECParticle.h`**

- **功能**：
  - 提供粒子信息的接口，包括粒子的动量分量、顶点位置、电荷等。
  - 提供工具函数，用于筛选特定粒子并提取其属性。
- **主要逻辑**：
  1. 定义粒子数据的结构（如 `RECParticle::All()`、`RECParticle::Extend()`等）。
  2. 提供` get_RECParticle_float_var` 和 `get_RECParticle_int_var` 函数，用于筛选粒子并提取其属性。
  3. `RECParticle::Extend()` 粒子数据的扩展结构（如添加 `theta`, `phi` 和 `p`）。

---

#### **3. `Math/RECParticleKinematic.h`和 `Math/RECParticleKinematic.cxx`**

- **功能**：
  - 提供粒子运动学变量的计算函数，包括 `theta`, `phi` 和 `p`

- **主要逻辑**：
  1. 定义 `RECParticletheta`：计算粒子的极角。
  2. 定义 `RECParticlephi`：计算粒子的方位角。
  3. 定义 `RECParticleP`：计算粒子的动量。

---

#### **4. core/FilesInPath.h`和 `core/FilesInPath.cxx`**

- **功能**：
  - 提供工具函数 `GetHipoFilesInPath`，用于递归搜索指定目录下的所有 `.hipo` 文件。

- **主要逻辑**：
  1. 使用 C++17 的 `std::filesystem` 遍历目录。
  2. 检查文件扩展名是否为 `.hipo`。
  3. 返回所有符合条件的文件路径。
- **依赖**：
  - C++17 标准库。

---

#### **5. `DrawHist/DrawAndSave.h` 和 `DrawHistDrawAndSave.cxx`**

- **功能**：
  - 提供绘制和保存直方图的工具函数。
- **主要逻辑**：
  1. 定义 `DrawAndSaveParticleHistograms` 函数，用于生成粒子的 `theta`、`phi`和 `p` 直方图。
  2. 动态生成直方图的名称和标题。
  3. 使用 `Draw1DHist` 和 `Save1DHist` 函数绘制并保存直方图为`.png`和`.root`。

---

## **使用方法**

#### **1. 编译项目**

在Jefferson Lab ifarm上依次输入以下进行编译：

```
module load clas12
```

```bash
git clone https://github.com/yijie086/DISANA
```

```
cd DISANA
```

```
mkdir build
```

```
cd build
```

```
cmake ..
```

```
make
```

这一步也可以:

```
make -j8
```

#### **2. 运行程序**

运行生成的可执行文件 `main101`，并指定包含 `.hipo` 文件的目录：

```bash
./main101 /path/to/hipo/files
```

#### **3. 输出结果**

- 筛选后的数据将保存到 `snapshot.root` 文件中。
- 生成的直方图将保存到 `test.root` 文件中。
- 还将看到输出直方图的`.png`文件。

---

### **注意事项**

1. 如果需要修改筛选条件，无需编辑 `EventCut`，修改方法参考`run101.C`
2. 如果需要添加新的粒子属性或运动学变量，可以扩展 `RECParticle` 和 `RECParticleKinematic`。
