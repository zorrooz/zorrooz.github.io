# gblog

<center>
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Vue 3](https://img.shields.io/badge/Vue-3-4fc08d)](https://vuejs.org/)
[![Vite](https://img.shields.io/badge/Vite-7-646cff)](https://vitejs.dev/)
[![GitHub stars](https://img.shields.io/github/stars/zorrooz/gblog?style=social)](https://github.com/zorrooz/gblog)
[![中文版本](https://img.shields.io/badge/中文版本-blue)](README.md)
</center>

<center>[English](README_EN.md) | [中文](README.md)</center>

gblog is a modern personal blog system built with Vue 3 + Vite, featuring bilingual support (Chinese/English), theme switching, Markdown article management, and more.

## 📋 Table of Contents

[TOC]

## ✨ Introduction

gblog is a fully-featured static blog system designed for personal knowledge management and content presentation. Built using Vue 3 Composition API and Vite as the build tool, it offers responsive design and a modern user experience.

### 🎯 Core Features

- **Bilingual Support**: Built-in Chinese and English language switching
- **Theme Switching**: Support for light/dark theme modes
- **Markdown Rendering**: Full Markdown support including code highlighting and math formulas
- **Content Organization**: Multi-dimensional content management with categories, resources, about page, etc.
- **Auto-translation**: Integrated AI translation tools for automatic content translation
- **Responsive Design**: Adapts to various screen sizes
- **Fast Loading**: Quick builds and hot updates powered by Vite

### 🛠️ Tech Stack

- **Frontend Framework**: Vue 3 (Composition API)
- **Build Tool**: Vite 7
- **State Management**: Pinia
- **Routing**: Vue Router 4
- **UI Framework**: Bootstrap 5
- **Theme Implementation**: Sass
- **Markdown Processing**: unified ecosystem (remark & rehype)
- **Internationalization**: vue-i18n
- **Deployment**: GitHub Pages

### 📁 Project Structure

```
src/
├── content-src/         # Source content files
│   ├── notes/           # Note articles
│   ├── projects/        # Project documentation
│   ├── topics/          # Research topics
│   ├── categories.yaml  # Category definitions
│   ├── about.yaml       # About page content
│   └── resources.yaml   # Resources page content
├── content/             # Generated JSON files
├── views/               # Page components
├── components/          # Reusable components
├── stores/              # Pinia state management
├── utils/               # Utility functions & generators
└── router/              # Routing configuration
```

## 🚀 Getting Started

### Prerequisites

- **Node.js**: `^20.19.0` or `>=22.12.0`

### Installation

1. Clone the repository

   ```
   git clone https://github.com/zorrooz/zorrooz.github.io.git
   cd zorrooz.github.io
   ```

2. Install dependencies

   ```
   npm install
   ```

3. Start the development server
   
   ```
   npm run dev
   ```

   The dev server will run at http://localhost:5173 by default.

4. Build for production

   ```
   npm run build
   ```

5. Preview the production build

   ```
   npm run preview
   ```

6. Deploy to GitHub Pages
   
   Update the URL in the `deploy` command in `package.json` to your GitHub Pages address, then run:

   ```
   npm run deploy
   ```

### 📖 Detailed Usage

gblog uses a pure Markdown + YAML metadata approach for content management. The system uses build-time generators to convert source files into optimized JSON files for fast content loading and navigation.

**Workflow:**

1. **Writing Phase**: Create Markdown files and YAML configurations in `src/content-src/`
2. **Build Phase**: The `prebuild` script automatically runs generators to convert content to JSON
3. **Runtime Phase**: Vue components load JSON files and render content

#### Category Content Management

The project supports three content types, defined in `src/content-src/categories.yaml`:

1. **Notes**

   For learning records and technical documentation:

   ```yaml
   - name: "category_identifier"
     title: "Display Title"
     desc: "Category description"
     date: "creation_date"
     categories:
       subcategory_key: "Subcategory Display Name"
       another_subcategory: "Subcategory Display Name"
       ...
   ```

2. **Projects**

   For code projects and practical cases:

   ```yaml
   - name: "project_identifier"
     title: "Project Display Name"
     desc: "Project description"
     github: "GitHub repository URL"
     date: "project_creation_date"
     categories:
       subcategory_key: "Subcategory Display Name"
       another_subcategory: "Subcategory Display Name"
       ...
   ```

3. **Topics**

   For research topics and academic content:
   
   ```yaml
   - name: "topic_identifier"
     title: "Topic Display Name"
     desc: "Research description"
     doi: "Academic identifier"
     date: "topic_end_date"
     categories:
       subcategory_key: "Subcategory Display Name"
       another_subcategory: "Subcategory Display Name"
       ...
   ```

#### Writing Markdown Files

This project uses pure Markdown with standard Markdown syntax support. Edit metadata in `src/content-src/categories.yaml`.

File organization structure:

```
src/content-src/
├── categories.yaml          # Category definitions
├── categories-en.yaml       # English category definitions
├── notes/                   # Notes directory
│   ├── category_identifier/
│   │   ├── subcategory_key/
│   │   │   └── article_name.md
│   │   └── another_subcategory/
│   │       └── article_name.md
│   └── another_category/
├── projects/                # Projects directory
│   ├── project_identifier/
│   │   ├── subcategory/
│   │   │   └── project_doc.md
│   │   └── README.md
│   └── other_project/
└── topics/                  # Topics directory
    ├── topic_identifier/
    │   ├── research_direction/
    │   │   └── research_paper.md
    │   └── experiment_data/
    └── other_topic/
```

Workflow:

1. Create or edit Markdown files in the `src/content-src` directory
2. Write your article content
3. Run `npm run dev` to see the results

#### Resources Page Configuration

The resources page uses a three-level hierarchy: Category → Subcategory → Resource Item, defined in `src/content-src/resources.yaml`:

```yaml
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

The about page uses a three-part structure: Personal Introduction + Content Blocks + Contact Information, defined in `src/content-src/about.yaml`:

```yaml
# Personal introduction (required)
introduction: "Your self-introduction here. Supports multi-line text."

# Content blocks (optional, multiple allowed)
section:
  - title: "Block Title"
    items:
      - item: "Project name or position"
        desc: "Detailed description"
      - item: "Another project"
        desc: "Additional information"

# Contact information (optional, multiple allowed)
contacts:
  - label: "Contact Type"
    value: "Display text"
    link: "Link address"
    icon: "Font Awesome icon class"
```

#### Using the Translation Tool

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
  npm run translate -- status src/content-src
  ```

- More usage options
  
  Use `npm run translate -- help` for detailed usage.

#### Content Generators

Content generators typically run automatically during build time. If you need to manually update all content, run:

```
npm run prebuild
```

The generator sequentially processes:

1. Generates notes, projects, and topics indexes
2. Category structure generation (depends on step 1)
3. Article list generation (depends on notes information)
4. Tag cloud generation (depends on article information)

## ⚙️ Acknowledgments

This project is personally developed by Zorrooz with assistance from multiple LLMs and AI IDEs. Special thanks to the following services:

- [deepseek](https://chat.deepseek.com/)
- [Qianwen](https://www.qianwen.com/)
- [Doubao](https://www.doubao.com/chat/)
- [Trae](https://www.trae.cn/)
- [CodeBuddy](https://www.codebuddy.ai/)
- [Visual Studio Code (GitHub Copilot)](https://code.visualstudio.com/)

## 📄 License

Copyright (c) 2025 Zorrooz. All rights reserved.

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
**If you find this project helpful, please consider giving it a Star ⭐. Thank you!**
