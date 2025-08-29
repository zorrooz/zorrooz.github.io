<template>
  <div class="markdown-content" v-html="renderedMarkdown"></div>
</template>

<script>
import md from 'markdown-it'
import hljs from 'highlight.js'
import 'highlight.js/styles/atom-one-dark.css' // 现代感的代码高亮主题
import 'highlight.js/lib/common' // 导入常见语言支持

// 自定义标题渲染规则
const renderTitle = (tokens, idx, options, env, self) => {
  const token = tokens[idx]
  const level = token.tag.slice(1) // 提取标题级别 (h1 -> 1)
  const titleText = tokens[idx + 1].content // 获取标题文本

  // 生成标题ID用于锚点链接 (将文本转为kebab-case)
  const titleId = titleText
    .toLowerCase()
    .replace(/[^a-z0-9\s]/g, '')
    .replace(/\s+/g, '-')

  // 添加锚点链接图标
  return `<${token.tag} id="${titleId}" class="md-title md-title-${level}">
            ${self.renderToken(tokens, idx, options)}
            <a href="#${titleId}" class="title-anchor" aria-label="Anchor for ${titleText}">
              #
            </a>
          </${token.tag}>`
}

export default {
  props: {
    content: {
      type: String,
      required: true,
      default: ''
    },
    options: {
      type: Object,
      default: () => ({
        breaks: true,
        html: false,
        linkify: true,
        typographer: true,
        highlight: (str, lang) => {
          // 代码高亮增强
          if (lang && hljs.getLanguage(lang)) {
            try {
              // 添加语言标签
              return `<div class="code-block">
                        <div class="code-language">${lang}</div>
                        <pre class="hljs"><code>${hljs.highlight(str, { language: lang, ignoreIllegals: true }).value}</code></pre>
                      </div>`
            } catch (err) {
              console.warn(`Highlight error for language ${lang}:`, err)
            }
          }

          // 自动检测语言并添加通用标签
          return `<div class="code-block">
                    <div class="code-language">code</div>
                    <pre class="hljs"><code>${hljs.highlightAuto(str).value}</code></pre>
                  </div>`
        }
      })
    }
  },
  computed: {
    renderedMarkdown() {
      if (typeof this.content !== 'string') {
        console.error('Markdown content must be a string')
        return '<div class="error">Invalid markdown content</div>'
      }

      const markdown = md(this.options)

      // 重写标题渲染规则 (h1-h6)
      for (let i = 1; i <= 6; i++) {
        markdown.renderer.rules[`heading_open_${i}`] = renderTitle
      }

      return markdown.render(this.content)
    }
  }
}
</script>

<style scoped>
.markdown-content {
  max-width: 800px;
  margin: 0 auto;
  padding: 2rem;
  line-height: 1.8;
  color: #333;
  font-family: 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
}

/* 标题样式增强 */
.md-title {
  position: relative;
  margin: 2rem 0 1rem;
  padding-bottom: 0.5rem;
  color: #2c3e50;
  font-weight: 600;
  transition: color 0.3s ease;
}

/* 不同级别标题的差异化样式 */
.md-title-1 {
  font-size: 2rem;
  border-bottom: 3px solid #3498db;
}

.md-title-2 {
  font-size: 1.75rem;
  border-bottom: 2px solid #3498db;
}

.md-title-3 {
  font-size: 1.5rem;
  border-bottom: 1px solid #3498db;
}

.md-title-4 {
  font-size: 1.25rem;
  color: #34495e;
}

.md-title-5 {
  font-size: 1.1rem;
  color: #34495e;
}

.md-title-6 {
  font-size: 1rem;
  color: #7f8c8d;
}

/* 标题锚点链接 */
.title-anchor {
  position: absolute;
  left: -1.5rem;
  top: 50%;
  transform: translateY(-50%);
  color: #bdc3c7;
  text-decoration: none;
  opacity: 0;
  transition: all 0.2s ease;
}

.md-title:hover .title-anchor {
  opacity: 1;
  color: #3498db;
}

.title-anchor:hover {
  text-decoration: none;
  color: #2980b9;
}

/* 代码块增强样式 */
.code-block {
  position: relative;
  margin: 1.5rem 0;
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
}

.code-language {
  padding: 0.4rem 1rem;
  background-color: #2d2d2d;
  color: #f0f0f0;
  font-size: 0.8rem;
  text-transform: uppercase;
  letter-spacing: 0.5px;
}

pre.hljs {
  margin: 0;
  padding: 1.2rem;
  background-color: #2d2d2d;
  border-radius: 0 0 8px 8px;
}

/* 行内代码样式 */
.markdown-content code:not(pre code) {
  padding: 0.2em 0.4em;
  margin: 0 0.2em;
  background-color: #f5f5f5;
  border-radius: 4px;
  font-family: 'Fira Code', 'SFMono-Regular', Consolas, monospace;
  font-size: 0.9em;
}

/* 其他基础样式保持不变 */
.markdown-content a {
  color: #3498db;
  text-decoration: none;
  transition: color 0.2s;
}

.markdown-content a:hover {
  color: #2980b9;
  text-decoration: underline;
}

.markdown-content ul,
.markdown-content ol {
  padding-left: 2rem;
  margin: 1rem 0;
}

.error {
  color: #e74c3c;
  padding: 1rem;
  border: 1px solid #fadbd8;
  border-radius: 4px;
  background-color: #fef5f5;
}
</style>
