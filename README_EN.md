# gblog

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/zorrooz/zorrooz.github.io)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)
[![GitHub stars](https://img.shields.io/github/stars/zorrooz/gblog?style=social)](https://github.com/zorrooz/zorrooz.github.io)
[![English Version](https://img.shields.io/badge/English-Version-blue)](README_EN.md)


[English](README_EN.md) | [中文](README.md)

gblog is a modern personal blog system built with Vue 3 + Vite, featuring bilingual support (English/Chinese), theme switching, Markdown article management, and more.

## ✨ Project Introduction

gblog is a fully-featured static blog system designed for personal knowledge management and content presentation. Developed with Vue 3's Composition API and Vite as the build tool, it supports responsive design and a modern user experience.

### 🎯 Core Features

- **Bilingual Support**: Built-in English and Chinese language switching
- **Theme Switching**: Supports light / dark / system-follow modes
- **Markdown Rendering**: Full Markdown support with code highlighting and mathematical formulas
- **Content Copying**: One-click copy for articles and tables (tables copied as TSV, paste-ready for Excel)
- **Content Organization**: Multi-dimensional content management with categories, resources, about page, etc.
- **Auto Translation**: Integrated AI translation tool for automatic English-Chinese content translation
- **Responsive Design**: Adapts to various screen sizes
- **Fast Loading**: Quick builds and hot updates powered by Vite

### 🛠️ Tech Stack

- **Frontend Framework**: Vue 3 (Composition API)
- **Build Tool**: Vite 7
- **State Management**: Pinia
- **Routing Management**: Vue Router 4
- **UI Framework**: Bootstrap 5
- **Theme Implementation**: Sass
- **Markdown Processing**: unified ecosystem (remark and rehype)
- **Internationalization**: vue-i18n
- **Deployment**: GitHub Pages

### 📁 Project Structure

Content lives on the `data` branch (local worktree at `../blog-data`); `main` branch holds code only:

```
src/  
├── views/               # Page components  
├── components/          # Reusable components  
├── stores/              # Pinia state management  
├── utils/               # Utility functions and generators  
└── router/              # Routing configuration  

# data branch (../blog-data)
content-src/             # Layer 1 src: hand-written Chinese sources (edit here only)
│   ├── notes/           # Note articles (only directory with markdown)  
│   ├── projects/        # Project metadata (yaml only, no articles)  
│   ├── topics/          # Topic metadata (yaml only, no articles)  
│   ├── categories.yaml  # Category definitions  
│   ├── about.yaml       # About page content  
│   └── resources.yaml   # Resources page content  
cache/                   # Layer 2 cache: machine-maintained persistent state (committed)
│   ├── en/              #   English translation layer (mirrors content-src, -en identity)
│   ├── tag-mapping.json #   zh→en tag mapping
│   └── .translate-state.json # incremental translation state
content/                 # Layer 3 final: generated JSON + html (regenerable, not committed)
```

## 🚀 Getting Started

### 🌍 Environment

- **Node.js**: `>=23.6.0` (Node native type stripping runs `.ts` generators directly)

### 💻 Installation

1. Clone the repository and attach the data worktree (sibling dir `../blog-data`)

   ```
   git clone https://github.com/zorrooz/zorrooz.github.io.git  
   cd zorrooz.github.io
   git worktree add ../blog-data data
   ```

2. Install dependencies

   ```
   npm install
   ```

3. Start the development server
   
   ```
   npm run dev
   ```

   The development server will run at http://localhost:5173 by default.

4. Build for production

   ```
   npm run build
   ```

5. Preview the production build

   ```
   npm run preview
   ```

6. Publish (after editing `../blog-data/content-src/**`)

   Commit and push the data branch; GitHub Actions rebuilds and deploys automatically:

   ```
   npm run deploy
   ```

### 📖 Detailed Usage

gblog uses a pure Markdown + YAML metadata approach for content management. The system converts source files into optimized JSON files through a build-time generator, enabling fast content loading and navigation.

**Workflow:**

1. **Writing Phase**: Create Markdown files and YAML configurations in `content-src/` on the data branch (`../blog-data/content-src/`)
2. **Build Phase**: The `prebuild` script automatically runs generators to convert content to JSON (`../blog-data/content/`)
3. **Runtime Phase**: Vue components load JSON files and render page content
4. **Publish Phase**: `npm run deploy` commits and pushes the data branch; CI rebuilds and deploys

#### Category Content Management

The project supports three content types, defined in `content-src/categories.yaml`:

1. **Notes (notes)**

   For learning records and technical documentation:

   ```
   - name: "category_identifier"  
     title: "Display Title"  
     desc: "Category description"  
     date: "creation_date"  
     categories:  
       subcategory_key: "Subcategory Display Name"  
       another_subcategory_key: "Subcategory Display Name"
       ...
   ```

2. **Projects (projects)**

   For code projects and practical cases:

   ```
   - name: "project_identifier"  
     title: "Project Display Name"  
     desc: "Project function description"  
     github: "GitHub repository URL"  
     date: "project_creation_date"  
     categories:  
       subcategory_key: "Subcategory Display Name"  
       another_subcategory_key: "Subcategory Display Name"
       ...
   ```

3. **Topics (topics)**

   For research topics and academic content:
   
   ```
   - name: "topic_identifier"  
     title: "Topic Display Name"  
     desc: "Research content description"  
     doi: "academic_identifier"  
     date: "topic_end_date"  
     categories:  
       subcategory_key: "Subcategory Display Name"  
       another_subcategory_key: "Subcategory Display Name"
       ...
   ```

#### Writing Markdown Files

This project uses pure Markdown with standard Markdown syntax support. Please edit metadata in `content-src/categories.yaml`.

File organization structure:

```
content-src/  
├── categories.yaml          # Category definitions  
├── notes/                   # Notes directory (only directory with markdown)  
│   ├── category_identifier/  
│   │   ├── subcategory/  
│   │   │   └── article_name/
│   │   │       └── article_name.md  
│   │   └── other_subcategory/  
│   │       └── article_name/
│   │           └── article_name.md  
│   └── other_category/  
```

> Note: **projects and topics are metadata-only (yaml)** — defined in `categories.yaml` (fields like github/doi/url/language/license), with no markdown articles; their cards link out to GitHub / DOI directly.

Workflow:

1. Create or edit Markdown files in `content-src` on the data branch (`../blog-data/content-src/`)
2. Write article content
3. Run `npm run dev` to view the results

#### Resources Page Configuration

The resources page uses a three-level hierarchy: Category → Subcategory → Resource Item, defined in `content-src/resources.yaml`:

```
# Top-level category  
- title: "Category Name"  
  children:  
    # Subcategory  
    - title: "Subcategory Name"  
      items:  
        # Specific resource  
        - name: "Resource Name"  
          url: "https://example.com"  
          desc: "Resource description"
```

#### About Page Configuration

The about page uses a three-section structure: Personal Introduction + Content Blocks + Contact Information, defined in `content-src/about.yaml`:

```
# Personal introduction (required)  
introduction: "Your self-introduction here, supports multi-line text"  
  
# Content blocks (optional, can have multiple)  
section:  
  - title: "Block Title"  
    items:  
      - item: "Project name or position"  
        desc: "Detailed description or explanation"  
      - item: "Another project"  
        desc: "Additional information"  
  
# Contact information (optional, can have multiple)  
contacts:  
  - label: "Contact Type"  
    value: "Display text"  
    link: "Link address"  
    icon: "Font Awesome icon class"
```

#### Using the Translation Tool

Before starting the translation, please complete the API Key configuration:

Create a new file `llmConfig.ts` in the `scripts` directory (it is gitignored and never committed) and add the following content:

```
export default {
  url: 'your_api_endpoint_url',  // Replace with the actual API endpoint URL
  apikey: 'your_api_key_here',   // Replace with the actual API key
  model: 'your_model_name',      // Replace with the model name to be used
};

```

After configuration is complete, you can use the following commands for translation operations:

- Incremental translation (recommended)

  ```
  npm run translate
  ```

- Translate only new files
  
  ```
  npm run translate -- translate --new
  ```

- Force re-translation of all files

  ```
  npm run translate -- translate --force
  ```

- Check translation status

  ```
  npm run translate
  ```

- More usage
  
  Use `npm run translate -- help` to view detailed usage.

After translation is complete, please manually trigger a build to update the content:

```
npm run prebuild
```

#### Content Generators

Normally, content generators run automatically during build time without manual intervention (`runAllGenerators`).
To manually update all content, run:

```
npm run prebuild
```

The generator will process in sequence:

1. Generate notes, projects, topics indexes
2. Category structure generation (depends on step 1)
3. Article list generation (depends on notes information)
4. Tag cloud generation (depends on article information)

## ⚙️ Acknowledgments

This project was developed personally by Zorrooz, with assistance from multiple LLMs and AI IDEs for coding. Thanks to the following services:

- [deepseek](https://chat.deepseek.com/)
- [千问](https://www.qianwen.com/)
- [豆包](https://www.doubao.com/chat/)
- [Trae](https://www.trae.cn/)
- [CodeBuddy](https://www.codebuddy.ai/)
- [Visual Studio Code (GitHub Copilot)](https://code.visualstudio.com/)


## 📄 License

Copyright (c) 2025-2026 Zorrooz. All rights reserved.

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
**If this project is helpful to you, please consider giving it a star ⭐. Thank you!**
