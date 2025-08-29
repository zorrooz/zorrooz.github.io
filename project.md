# 项目设计

我要创建一个blog项目。网站风格为现代化的简约风格，主题颜色为#047AFF。项目结构：主页+归档页+项目页+标签页+资源页。所有页面分为三个部分：header footer 主体。一个核心：md网页转换器。header包含：最左侧依次为logo（点击可以跳转主页） 归档按钮 项目按钮 标签 资源按钮，右边从左到右依次为搜索图标（点击后弹出搜索框在主体部分最上方） 语言切换图标（中文和英文） 主题切换图标（深色和浅色）。footer包含网站的创建信息和本人github+email。网站核心是md显示：包含正常的md显示，右侧为on this page的目录。对于一篇文章，包含2个部分：js配置（包含创建时间、归类、封面、标签等元信息）和md原始文件，md第一段默认为内容预览。主页布局：左（上方卡片：头像、本人简介、[文章数量、分类、标签、字数]的统计信息，下方卡片：文章一级分类），中（一条一条文章卡片，下方可以点击下一页等），右（上方卡片：最近文章的条目，包含日期和2行内容预览。下方卡片：依据计数的热门标签列表，不要显示太多）。归档页面：依旧保持主页的左中右布局，只是中间改为自定义的树状文章分类列表，每一个一级分类的文章都会形成以下结构：book chapter article，book是一个一级分类，chapter是二级分类，article是三级分类，因此点击任意文章进去后，左侧显示整个book及其章节内容，右侧显示on this page，类似于开源项目文档的布局。项目页面：依旧保持主页的左中右布局，只是中间改为自定义的项目卡片式列表，大体分为两类项目：干实验项目+湿实验项目，采用卡片式布局，一个项目一个卡片横向布局，每行最多三个卡片，卡片上面是项目的基本信息，分为两类：干实验和湿实验，不要直接写干实验和湿实验，你换个更合适的名字。每个项目点进去以后和归档点进去布局差不多，只不过左侧目录为项目文档目录，右侧依然为on this page。标签页面：依旧保持主页的左中右布局，只是中间改为标签排布：按照拼音排序，其大小是标签下的文章计数，点击一个标签即可列出该标签下的全部文章。资源页面：依旧保持主页的左中右布局，只是中间改为个性化的资源快速导航链接。

```
blog-project/
├── public/                  # 静态资源
│   └── favicon.ico
├── src/
│   ├── router/              # 路由配置中心
│   │   └── index.js         # 路由定义和导航守卫
│   ├── assets/              # 静态资源
│   │   ├── css/
│   │   │   ├── main.css     # 全局样式
│   │   │   └── theme.css    # 深/浅色主题变量
│   │   └── images/          # 图片资源
│   ├── components/
│   │   ├── layout/
│   │   │   ├── AppHeader.vue # 头部组件（含导航/搜索/切换）
│   │   │   └── AppFooter.vue # 底部组件
│   │   ├── content/
│   │   │   ├── MarkdownRenderer.vue # MD解析核心组件
│   │   │   ├── TOC.vue      # "On This Page"目录组件
│   │   │   └── ArticleCard.vue
│   │   └── widgets/
│   │       ├── CategoryTree.vue # 归档树形结构
│   │       ├── TagCloud.vue  # 标签云
│   │       └── ProjectGrid.vue # 项目卡片网格
│   ├── content/             # MD内容中心
│   │   ├── posts/           # 文章（含元数据）
│   │   │   └── 2023-10-01-example.md
│   │   ├── projects/        # 项目文档
│   │   └── resources.md     # 资源页内容
│   ├── stores/              # Pinia状态管理
│   │   ├── themeStore.js    # 主题切换
│   │   ├── langStore.js     # 语言切换
│   │   └── searchStore.js   # 搜索功能
│   ├── utils/
│   │   ├── markdownParser.js # MD解析器
│   │   ├── contentLoader.js # 内容加载工具
│   │   └── pinyinSorter.js  # 中文标签排序
│   ├── views/
│   │   ├── Home.vue         # 主页（三栏布局）
│   │   ├── Archive.vue      # 归档页
│   │   ├── Projects.vue     # 项目页（干/湿实验分类）
│   │   ├── Tags.vue         # 标签页
│   │   ├── Resources.vue    # 资源页
│   │   └── PostDetail.vue   # 文章/项目详情页
│   ├── App.vue
│   └── main.js
├── vite.config.js           # Vite配置（含MD插件）
└── package.json
```
