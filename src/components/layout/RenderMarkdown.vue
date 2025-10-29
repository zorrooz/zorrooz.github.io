<template>
  <div class="markdown-body" v-html="renderedMarkdown" ref="markdownContainer"></div>
</template>

<script>
import { renderMarkdown } from '@/utils/markdownProcessor'

const assetModules = import.meta.glob('../../content-src/**/*.{png,jpg,jpeg,gif,svg,webp}', { as: 'url', eager: true });

export default {
  props: {
    rawMarkdown: { type: String, default: '' },
    articlePath: { type: String, default: '' },
    articleTitle: { type: String, default: '' }
  },
  data() { return { renderedMarkdown: '' } },
  watch: {
    rawMarkdown: { handler: 'initRender', immediate: true, async: true }
  },
  methods: {
    async initRender(markdownContent) {
      const processedMarkdown = this.rewriteImageLinks(markdownContent, this.articlePath)
      this.renderedMarkdown = await renderMarkdown(processedMarkdown)
      this.$nextTick(() => {
        this.$emit('markdown-rendered')
        this.enhanceCodeBlocks()
        this.enhanceHeadings()
      })
    },

    rewriteImageLinks(md, articlePath) {
      try {
        const articleDir = articlePath.replace(/^[./]*/, '').replace(/\.md$/, '').split('/').slice(0, -1).join('/');
        
        const toAssetUrl = (relPath) => {
          if (/^(https?:)?\/\//i.test(relPath) || relPath.startsWith('/')) return relPath;
          const parts = (articleDir + '/' + relPath).split('/').filter(p => p && p !== '.');
          const stack = [];
          parts.forEach(p => p === '..' ? stack.pop() : stack.push(p));
          const normalized = stack.join('/');
          // 尝试多种可能的路径格式
          const candidateKeys = [
            `../../content-src/${normalized}`,
            `../content-src/${normalized}`,
            `content-src/${normalized}`
          ];
          for (const key of candidateKeys) {
            if (assetModules[key]) return assetModules[key];
          }
          return relPath;
        };

        return md
          .replace(/!\[([^\]]*)\]\(([^)]+)\)/g, (m, alt, src) => {
            const clean = src.trim().replace(/^<|>|&/g, '');
            return `![${alt}](${toAssetUrl(clean)})`;
          })
          .replace(/<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi, (m, pre, src, post) => 
            `<img ${pre}src="${toAssetUrl(src.trim())}"${post}>`);
      } catch (e) {
        console.warn('rewriteImageLinks failed', e);
        return md;
      }
    },

    enhanceHeadings() {
      const container = this.$refs.markdownContainer;
      if (!container) return;
      this.cleanDuplicateH1(container);
      this.addAnchorLinks(container);
    },

    cleanDuplicateH1(container) {
      if (!this.articleTitle) return;
      const pageTitle = this.articleTitle.trim().toLowerCase();
      container.querySelectorAll('h1').forEach(h1 => {
        const h1Text = h1.textContent.trim().toLowerCase();
        if (h1Text === pageTitle) h1.remove();
        else h1.replaceWith(Object.assign(document.createElement('h2'), {
          ...Object.fromEntries(Array.from(h1.attributes).map(attr => [attr.name, attr.value])),
          innerHTML: h1.innerHTML
        }));
      });
    },

    addAnchorLinks(container) {
      const scrollToHeading = (heading, anchorBtn) => {
        const targetTop = window.scrollY + heading.getBoundingClientRect().top - 8;
        window.scrollTo({ top: Math.max(0, targetTop), behavior: 'smooth' });
        setTimeout(() => anchorBtn.blur(), 300);
      };

      container.querySelectorAll('h2, h3, h4, h5, h6').forEach(heading => {
        heading.querySelector('.heading-anchor')?.remove();
        const anchorBtn = Object.assign(document.createElement('button'), {
          type: 'button',
          className: 'heading-anchor',
          textContent: '#',
          ariaLabel: '置顶当前标题',
          tabIndex: 0,
          ariaHidden: 'false'
        });

        anchorBtn.addEventListener('click', (e) => { 
          e.stopPropagation(); 
          scrollToHeading(heading, anchorBtn); 
        });
        anchorBtn.addEventListener('keydown', (e) => {
          if (e.key === 'Enter' || e.key === ' ') { 
            e.preventDefault(); 
            scrollToHeading(heading, anchorBtn); 
          }
        });

        heading.appendChild(anchorBtn);
      });
    },

    enhanceCodeBlocks() {
      const container = this.$refs.markdownContainer
      if (!container) return

      container.querySelectorAll('pre').forEach(pre => {
        if (pre.querySelector('.code-block-header')) return
        const code = pre.querySelector('code')
        if (!code) return

        const language = (code.className.match(/language-(\w+)/) || [, 'text'])[1]
        const header = document.createElement('div')
        header.className = 'code-block-header'

        // 语言标签
        const langLabel = document.createElement('span')
        langLabel.className = 'code-language'
        langLabel.textContent = language

        // 复制按钮
        const copyButton = document.createElement('button')
        copyButton.type = 'button'
        copyButton.className = 'copy-button btn-icon'
        copyButton.setAttribute('aria-label', '复制代码')
        copyButton.innerHTML = `
          <svg width="18" height="18" viewBox="0 0 14 14" fill="currentColor">
            <path d="M3 2C2.44772 2 2 2.44772 2 3V9C2 9.55228 2.44772 10 3 10H9C9.55228 10 10 9.55228 10 9V3C10 2.44772 9.55228 2 9 2H3ZM1 3C1 1.89543 1.89543 1 3 1H9C10.1046 1 11 1.89543 11 3V9C11 10.1046 10.1046 11 9 11H3C1.89543 11 1 10.1046 1 9V3Z"/>
            <path d="M5 4C4.44772 4 4 4.44772 4 5V11C4 11.5523 4.44772 12 5 12H11C11.5523 12 12 11.5523 12 11V5C12 4.44772 11.5523 4 11 4H5Z"/>
          </svg>
        `
        copyButton.addEventListener('click', () => this.copyToClipboard(code.textContent, copyButton))

        // 组装DOM
        header.append(langLabel, copyButton)
        const wrapper = document.createElement('div')
        wrapper.className = 'code-block-wrapper'
        pre.parentNode.insertBefore(wrapper, pre)
        wrapper.append(header, pre)
      })
    },

    async copyToClipboard(text, button) {
      try {
        await navigator.clipboard.writeText(text)
      } catch (err) {
        console.error('复制失败:', err)
        // 降级方案
        const textArea = document.createElement('textarea')
        textArea.value = text
        document.body.appendChild(textArea)
        textArea.select()
        document.execCommand('copy')
        document.body.removeChild(textArea)
      } finally {
        this.showCopyFeedback(button)
      }
    },

    showCopyFeedback(button) {
      const originalColor = button.style.color
      button.style.color = 'var(--app-success)'
      setTimeout(() => button.style.color = originalColor, 1000)
    }
  }
}
</script>

<style>
/* 代码字体定义 */
@font-face {
  font-family: 'CodeFont';
  src: local('Agave Regular'), local('Agave-Regular'), url('@/assets/fonts/Agave-Regular.ttf') format('truetype');
  font-weight: normal;
  font-style: normal;
  unicode-range: U+0000-00FF, U+0131, U+0152-0153, U+02C6, U+02DA, U+02DC, U+2000-206F, U+2074, U+20AC, U+2212, U+2215;
}
@font-face {
  font-family: 'CodeFont';
  src: local('Source Han Sans SC Regular'), local('SourceHanSansSC-Regular'), url('@/assets/fonts/SourceHanSansSC-Regular.otf') format('opentype');
  font-weight: normal;
  font-style: normal;
  unicode-range: U+4E00-9FFF, U+3400-4DBF, U+F900-FAFF, U+20000-2A6DF, U+2A700-2B73F, U+2B740-2B81F, U+2B820-2CEAF, U+2CEB0-2EBEF, U+30000-3134F;
}

/* Markdown基础样式 */
.markdown-body {
  box-sizing: border-box;
  width: 100%;
  padding: 0;
  color: var(--app-text);
  line-height: 1.8;
}

/* 列表样式 */
.markdown-body ul, .markdown-body ol { padding-left: 1.5em; margin-bottom: 1em; }
.markdown-body ul { list-style-type: disc; }
.markdown-body ol { list-style-type: decimal; }
.markdown-body li { margin-bottom: 0.25em; padding-left: 0; }
.markdown-body ul li::marker, .markdown-body ol li::marker { color: var(--app-primary-mid); }

/* 链接样式 */
.markdown-body a { color: var(--app-primary); text-decoration: none; transition: color 0.2s ease; }
.markdown-body a:hover { color: var(--app-primary); text-decoration: underline; }
.markdown-body a[href^="http"]:not([href*="localhost"]):not([href*="127.0.0.1"]) {
  color: var(--app-markdown-external-link-color);
  font-weight: 700;
  text-decoration: none;
}
.markdown-body a[href^="http"]:not([href*="localhost"]):not([href*="127.0.0.1"]):hover { text-decoration: underline; }

/* 行内代码样式 */
.markdown-body code:not(pre code) {
  background-color: var(--app-markdown-code-bg);
  color: var(--app-markdown-code-color);
  padding: 0.2em 0.4em;
  border-radius: 0.25em;
  font-size: 0.9em;
  font-family: 'CodeFont', 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
}

/* 代码块样式 */
.markdown-body pre {
  background-color: var(--app-markdown-code-bg);
  border-radius: 0 0 0.5em 0.5em;
  padding: 1em;
  margin: 0;
  overflow-x: auto;
  font-family: 'CodeFont', 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
}
.markdown-body pre code { background: none; padding: 0; color: inherit; font-family: inherit; }

/* 代码块包装器 */
.markdown-body .code-block-wrapper {
  position: relative;
  margin: 1em 0;
  border-radius: 0.5em;
  overflow: hidden;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
}

/* 代码块头部 */
.markdown-body .code-block-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0.4em 1em;
  background-color: var(--app-markdown-code-header-bg, rgba(0, 0, 0, 0.05));
  border-bottom: 1px solid var(--app-border);
}

/* 语言标签 */
.markdown-body .code-language {
  font-size: 0.85em;
  font-weight: 500;
  color: var(--app-code-header-text);
  font-family: 'Inter', 'Roboto', 'Helvetica Neue', Arial, sans-serif;
  text-transform: lowercase;
  letter-spacing: 0.5px;
  padding: 0;
}

/* 复制按钮 */
.markdown-body .copy-button {
  background: transparent;
  border: none;
  padding: 0.25em;
  cursor: pointer;
  color: var(--app-text-muted);
  display: flex;
  align-items: center;
  justify-content: center;
}
.markdown-body .copy-button:hover { color: var(--app-primary); }

/* 语法高亮 */
.markdown-body .hljs { background: var(--app-markdown-code-bg); font-family: 'CodeFont', 'Monaco', 'Menlo', 'Ubuntu Mono', monospace; }

/* 图片样式 */
.markdown-body img {
  max-width: 100%;
  height: auto;
  display: block;
  margin: 1em auto;
  border-radius: 0.25em;
}

/* 段落和标题 */
.markdown-body p { margin-top: 0; margin-bottom: 1em; }
.markdown-body h1, .markdown-body h2, .markdown-body h3, .markdown-body h4, .markdown-body h5, .markdown-body h6 {
  font-weight: 600;
  line-height: 1.25;
  color: var(--app-text-emphasis);
}

.markdown-body h1 {
  margin-top: 1.2em;
  margin-bottom: 0.5em;
}

.markdown-body h2 {
  margin-top: 1em;
  margin-bottom: 0.4em;
}

.markdown-body h3 {
  margin-top: 0.8em;
  margin-bottom: 0.3em;
}

.markdown-body h4 {
  margin-top: 0.6em;
  margin-bottom: 0.2em;
}

.markdown-body h5 {
  margin-top: 0.5em;
  margin-bottom: 0.15em;
}

.markdown-body h6 {
  margin-top: 0.4em;
  margin-bottom: 0.1em;
}
.markdown-body h1 { font-size: 2em; }
.markdown-body h2 { font-size: 1.5em; }
.markdown-body h3 { font-size: 1.25em; }
.markdown-body h4 { font-size: 1.1em; }
.markdown-body h5 { font-size: 1em; }
.markdown-body h6 { font-size: 0.9em; color: var(--app-text-emphasis); }

/* 标题锚点链接 */
.markdown-body h2, .markdown-body h3, .markdown-body h4, .markdown-body h5, .markdown-body h6 {
  position: relative; 
  padding-right: 0; 
  display: inline-block; 
  width: 100%;
}

.markdown-body .heading-anchor {
  position: relative;
  display: inline-block;
  margin-left: 0.3em;
  color: var(--app-text-muted);
  font-size: 0.9em;
  font-weight: 400;
  opacity: 0.6;
  transition: opacity 0.2s ease, color 0.2s ease;
  cursor: pointer;
  background: none;
  border: none;
  border-radius: 0;
  box-shadow: none;
  padding: 0;
  line-height: 1;
}

.markdown-body h2:hover .heading-anchor, 
.markdown-body h3:hover .heading-anchor, 
.markdown-body h4:hover .heading-anchor, 
.markdown-body h5:hover .heading-anchor, 
.markdown-body h6:hover .heading-anchor, 
.heading-anchor:focus {
  opacity: 1;
  color: var(--app-primary-mid);
}

.markdown-body .heading-anchor:focus {
  outline: none;
  outline-offset: 0;
}

/* 引用块 */
.markdown-body blockquote {
  border-left: 4px solid var(--app-border);
  padding-left: 1em;
  margin: 1em 0;
  color: var(--app-text-muted);
  font-style: normal;
}

/* 表格 */
.markdown-body table { border-collapse: collapse; width: 100%; margin: 1em 0; }
.markdown-body th, .markdown-body td { border: 1px solid var(--app-border); padding: 0.5em; text-align: left; }
.markdown-body th { background-color: var(--app-bg-light); font-weight: 600; }
</style>
