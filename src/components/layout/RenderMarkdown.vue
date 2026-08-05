<!-- RenderMarkdown.vue -->
<template>
  <div class="markdown-body" v-html="renderedMarkdown" ref="markdownContainer"></div>
</template>

<script setup lang="ts">
import { nextTick, ref, useTemplateRef, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { copyText } from '@/utils/clipboard'

const { t } = useI18n()

// @data 别名指向数据分支；文章图片均从数据目录解析（zh 在 content-src，镜像层 cache/en）
const assetModules = import.meta.glob(
  '@data/{content-src,cache/en}/**/*.{png,jpg,jpeg,gif,svg,webp}',
  { query: '?url', import: 'default', eager: true },
)

const COPY_ICON_SVG = `
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M3 2C2.44772 2 2 2.44772 2 3V9C2 9.55228 2.44772 10 3 10H9C9.55228 10 10 9.55228 10 9V3C10 2.44772 9.55228 2 9 2H3ZM1 3C1 1.89543 1.89543 1 3 1H9C10.1046 1 11 1.89543 11 3V9C11 10.1046 10.1046 11 9 11H3C1.89543 11 1 10.1046 1 9V3Z"/>
    <path d="M5 4C4.44772 4 4 4.44772 4 5V11C4 11.5523 4.44772 12 5 12H11C11.5523 12 12 11.5523 12 11V5C12 4.44772 11.5523 4 11 4H5Z"/>
  </svg>
`

const CHECK_ICON_SVG = `
  <svg width="16" height="16" viewBox="0 0 14 14" fill="currentColor">
    <path d="M11.3536 3.64645C11.5488 3.84171 11.5488 4.15829 11.3536 4.35355L5.35355 10.3536C5.15829 10.5488 4.84171 10.5488 4.64645 10.3536L2.64645 8.35355C2.45118 8.15829 2.45118 7.84171 2.64645 7.64645C2.84171 7.45118 3.15829 7.45118 3.35355 7.64645L5 9.29289L10.6464 3.64645C10.8417 3.45118 11.1583 3.45118 11.3536 3.64645Z"/>
  </svg>
`

const props = withDefaults(
  defineProps<{
    rawMarkdown?: string
    articlePath?: string
    articleTitle?: string
  }>(),
  {
    rawMarkdown: '',
    articlePath: '',
    articleTitle: '',
  },
)

const emit = defineEmits(['markdown-rendered'])

const renderedMarkdown = ref('')
const markdownContainer = useTemplateRef<HTMLElement>('markdownContainer')

async function initRender(htmlContent: string) {
  const processedHtml = rewriteImageLinks(htmlContent, props.articlePath)
  renderedMarkdown.value = processedHtml
  await nextTick()
  emit('markdown-rendered')
  enhanceCodeBlocks()
  enhanceTables()
  enhanceHeadings()
}

function rewriteImageLinks(html: string, articlePath: string) {
  try {
    const articleDir = articlePath
      .replace(/^[./]*/, '')
      .replace(/\.md$/, '')
      .split('/')
      .slice(0, -1)
      .join('/')

    const toAssetUrl = (relPath: string) => {
      if (/^(https?:)?\/\//i.test(relPath) || relPath.startsWith('/')) return relPath
      const parts = (articleDir + '/' + relPath).split('/').filter((p) => p && p !== '.')
      const stack: string[] = []
      parts.forEach((p) => (p === '..' ? stack.pop() : stack.push(p)))
      const normalized = stack.join('/')
      const candidateKeys = [`@data/content-src/${normalized}`, `${normalized}`]
      for (const key of candidateKeys) {
        const matched = Object.keys(assetModules).find(
          (k) => k.endsWith(`/${normalized}`) || k === key,
        )
        if (matched) return assetModules[matched]
      }
      return relPath
    }

    return html.replace(
      /<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi,
      (_m, pre, src, post) => `<img ${pre}src="${toAssetUrl(src.trim())}"${post}>`,
    )
  } catch (e) {
    console.warn('rewriteImageLinks failed', e)
    return html
  }
}

function enhanceHeadings() {
  const container = markdownContainer.value
  if (!container) return
  cleanDuplicateH1(container)
  addAnchorLinks(container)
}

function cleanDuplicateH1(container: HTMLElement) {
  if (!props.articleTitle) return
  const pageTitle = props.articleTitle.trim().toLowerCase()
  container.querySelectorAll('h1').forEach((h1) => {
    const h1Text = h1.textContent.trim().toLowerCase()
    if (h1Text === pageTitle) h1.remove()
    else
      h1.replaceWith(
        Object.assign(document.createElement('h2'), {
          ...Object.fromEntries(Array.from(h1.attributes).map((attr) => [attr.name, attr.value])),
          innerHTML: h1.innerHTML,
        }),
      )
  })
}

function addAnchorLinks(container: HTMLElement) {
  const scrollToHeading = (heading: Element, anchorBtn: HTMLButtonElement) => {
    const targetTop = window.scrollY + heading.getBoundingClientRect().top - 88
    window.scrollTo({ top: Math.max(0, targetTop), behavior: 'smooth' })
    setTimeout(() => anchorBtn.blur(), 300)
  }

  container.querySelectorAll('h2, h3, h4, h5, h6').forEach((heading) => {
    heading.querySelector('.heading-anchor')?.remove()
    const anchorBtn = Object.assign(document.createElement('button'), {
      type: 'button',
      className: 'heading-anchor',
      textContent: '#',
      ariaLabel: t('anchorHeading'),
      tabIndex: 0,
      ariaHidden: 'false',
    })

    anchorBtn.addEventListener('click', (e) => {
      e.stopPropagation()
      scrollToHeading(heading, anchorBtn)
    })
    anchorBtn.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault()
        scrollToHeading(heading, anchorBtn)
      }
    })

    heading.appendChild(anchorBtn)
  })
}

function enhanceCodeBlocks() {
  const container = markdownContainer.value
  if (!container) return

  container.querySelectorAll('pre').forEach((pre) => {
    if (pre.querySelector('.code-block-header')) return
    const code = pre.querySelector('code')
    if (!code) return

    const language = (code.className.match(/language-(\w+)/) || ['', 'text'])[1]
    const header = document.createElement('div')
    header.className = 'code-block-header d-flex align-items-center justify-content-between'

    const langLabel = document.createElement('span')
    langLabel.className = 'code-language'
    langLabel.textContent = language

    const copyButton = document.createElement('button')
    copyButton.type = 'button'
    copyButton.className = 'copy-button btn-icon d-flex align-items-center justify-content-center'
    copyButton.setAttribute('aria-label', t('copyCode'))
    copyButton.innerHTML = COPY_ICON_SVG
    copyButton.addEventListener('click', () => copyToClipboard(code.textContent ?? '', copyButton))

    header.append(langLabel, copyButton)
    const wrapper = document.createElement('div')
    wrapper.className = 'code-block-wrapper'
    pre.parentNode?.insertBefore(wrapper, pre)
    wrapper.append(header, pre)
  })
}

function enhanceTables() {
  const container = markdownContainer.value
  if (!container) return

  container.querySelectorAll('table').forEach((table) => {
    if (table.closest('.table-copyable')) return

    const wrapper = document.createElement('div')
    wrapper.className = 'table-copyable'

    const copyButton = document.createElement('button')
    copyButton.type = 'button'
    copyButton.className = 'table-copy-btn'
    copyButton.setAttribute('aria-label', t('copyTable'))
    copyButton.innerHTML = COPY_ICON_SVG
    copyButton.addEventListener('click', () => copyTableToClipboard(table, copyButton))

    wrapper.append(copyButton)
    table.parentNode?.insertBefore(wrapper, table)
    wrapper.append(table)
  })
}

function copyTableToClipboard(table: HTMLTableElement, button: HTMLButtonElement) {
  const rows = Array.from(table.querySelectorAll('tr'))
  const lines = rows.map((tr) =>
    Array.from(tr.querySelectorAll('th, td'))
      .map((cell) => (cell.textContent || '').trim().replace(/\s+/g, ' '))
      .join('\t'),
  )
  copyToClipboard(lines.join('\n'), button)
}

async function copyToClipboard(text: string, button: HTMLButtonElement) {
  try {
    await copyText(text)
  } catch (err) {
    console.error(t('copyFailed'), err)
  } finally {
    showCopyFeedback(button)
  }
}

function showCopyFeedback(button: HTMLButtonElement) {
  const original = button.innerHTML
  button.style.color = 'var(--primary)'
  button.innerHTML = CHECK_ICON_SVG
  setTimeout(() => {
    button.innerHTML = original
    button.style.color = ''
  }, 1200)
}

watch(
  () => props.rawMarkdown,
  (value: string) => initRender(value),
  { immediate: true },
)
</script>

<style>
@font-face {
  font-family: 'CodeFont';
  src:
    local('Agave Regular'),
    local('Agave-Regular'),
    url('@/assets/fonts/Agave-Regular.ttf') format('truetype');
  font-weight: normal;
  font-style: normal;
  unicode-range:
    U+0000-00FF, U+0131, U+0152-0153, U+02C6, U+02DA, U+02DC, U+2000-206F, U+2074, U+20AC, U+2212,
    U+2215;
}

@font-face {
  font-family: 'CodeFont';
  src:
    local('Source Han Sans SC Regular'),
    local('SourceHanSansSC-Regular'),
    url('@/assets/fonts/SourceHanSansSC-Regular.otf') format('opentype');
  font-weight: normal;
  font-style: normal;
  unicode-range:
    U+4E00-9FFF, U+3400-4DBF, U+F900-FAFF, U+20000-2A6DF, U+2A700-2B73F, U+2B740-2B81F,
    U+2B820-2CEAF, U+2CEB0-2EBEF, U+30000-3134F;
}

.markdown-body {
  box-sizing: border-box;
  width: 100%;
  padding: 0;
  color: var(--app-text);
  line-height: 1.75;
}

.markdown-body ul,
.markdown-body ol {
  padding-left: 1.5em;
  margin-bottom: 1em;
}

.markdown-body ul {
  list-style-type: disc;
}

.markdown-body ol {
  list-style-type: decimal;
}

.markdown-body li {
  margin-bottom: 0.25em;
  padding-left: 0;
}

.markdown-body ul li::marker,
.markdown-body ol li::marker {
  color: var(--app-primary-mid);
}

.markdown-body a {
  color: var(--app-link);
  text-decoration: none;
  transition: color 0.14s ease;
}

.markdown-body a:hover {
  color: var(--app-link);
  text-decoration: underline;
}

.markdown-body a[href^='http']:not([href*='localhost']):not([href*='127.0.0.1']) {
  color: var(--app-markdown-external-link-color);
  font-weight: 600;
  text-decoration: none;
}

.markdown-body a[href^='http']:not([href*='localhost']):not([href*='127.0.0.1'])::after {
  content: '\f08e';
  font-family: 'Font Awesome 6 Free';
  font-weight: 900;
  font-size: 0.72em;
  margin-left: 0.35em;
  opacity: 0.5;
  transition: opacity 0.14s ease;
}

.markdown-body a[href^='http']:not([href*='localhost']):not([href*='127.0.0.1']):hover::after {
  opacity: 1;
}

.markdown-body a[href^='http']:not([href*='localhost']):not([href*='127.0.0.1']):hover {
  text-decoration: underline;
}

.markdown-body kbd {
  font-family: var(--font-mono);
  font-size: 0.82em;
  font-weight: 600;
  color: var(--fg-2);
  background: var(--tint);
  border: 1px solid var(--line);
  border-bottom-width: 2px;
  border-radius: 4px;
  padding: 0.1em 0.4em;
}

.markdown-body hr {
  border: none;
  height: 1px;
  background: var(--line);
  margin: 2em 0;
}

.markdown-body li.task-list-item {
  list-style: none;
  margin-left: -1.4em;
}

.markdown-body .task-list-item-checkbox {
  accent-color: var(--seed-primary);
  margin-right: 0.5em;
}

.markdown-body code:not(pre code) {
  background-color: var(--app-markdown-code-bg);
  color: var(--app-markdown-code-color);
  padding: 0.18em 0.45em;
  border-radius: var(--radius-sm);
  font-size: 0.88em;
  border: 1px solid var(--line);
  font-family: 'CodeFont', 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
}

.markdown-body pre {
  background-color: var(--app-markdown-code-bg);
  padding: 1.1em 1.25em;
  margin: 0;
  overflow-x: auto;
  font-size: 13.5px;
  line-height: 1.65;
  font-family: 'CodeFont', 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
}

.markdown-body pre code {
  background: none;
  padding: 0;
  color: inherit;
  font-family: inherit;
}

.markdown-body .code-block-wrapper {
  position: relative;
  margin: 1.5em 0;
  border-radius: var(--radius);
  overflow: hidden;
  border: 1px solid var(--line);
  box-shadow: var(--shadow-soft);
}

.markdown-body .table-copyable {
  position: relative;
}

.markdown-body .table-copy-btn {
  position: absolute;
  top: 8px;
  right: 8px;
  z-index: 1;
  display: flex;
  align-items: center;
  justify-content: center;
  width: 30px;
  height: 30px;
  border: 1px solid var(--line);
  border-radius: var(--radius-btn);
  background: var(--surface);
  color: var(--fg-3);
  cursor: pointer;
  opacity: 0;
  transition:
    opacity 0.14s ease,
    color 0.14s ease,
    background-color 0.14s ease;
}

@media (hover: hover) {
  .markdown-body .table-copyable:hover .table-copy-btn,
  .markdown-body .table-copy-btn:focus-visible {
    opacity: 1;
  }
}

@media (hover: none) {
  .markdown-body .table-copy-btn {
    opacity: 0.8;
  }
}

.markdown-body .table-copy-btn:hover,
.markdown-body .table-copy-btn:focus-visible {
  color: var(--primary);
  background: var(--tint);
  border-color: color-mix(in srgb, var(--primary) 30%, transparent);
}

.markdown-body .code-block-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0.55em 1.1em;
  background-color: var(--app-markdown-code-header-bg);
  border-bottom: 1px solid var(--line);
}

.markdown-body .code-language {
  font-size: 0.72rem;
  font-weight: 600;
  color: var(--primary);
  font-family: var(--font-mono);
  text-transform: lowercase;
  letter-spacing: 0.6px;
  padding: 0;
}

.markdown-body .copy-button {
  background: transparent;
  border: none;
  padding: 0.2em;
  cursor: pointer;
  color: var(--app-text-muted);
}

.markdown-body .copy-button:hover {
  color: var(--app-primary);
}

.markdown-body pre code.hljs {
  background: var(--app-markdown-code-bg);
  padding: 0;
  font-family: 'CodeFont', 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
}

.markdown-body p {
  margin-top: 0;
  margin-bottom: 1em;
}

.markdown-body h1,
.markdown-body h2,
.markdown-body h3,
.markdown-body h4,
.markdown-body h5,
.markdown-body h6 {
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

.markdown-body h1 {
  font-size: 2em;
}

.markdown-body h2 {
  font-size: 1.5em;
}

.markdown-body h3 {
  font-size: 1.25em;
}

.markdown-body h4 {
  font-size: 1.1em;
}

.markdown-body h5 {
  font-size: 1em;
}

.markdown-body h6 {
  font-size: 0.9em;
  color: var(--app-text-emphasis);
}

.markdown-body h2,
.markdown-body h3,
.markdown-body h4,
.markdown-body h5,
.markdown-body h6 {
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
  font-size: 0.85em;
  font-weight: 400;
  opacity: 0.4;
  transition:
    opacity 0.14s ease,
    color 0.14s ease;
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

.markdown-body blockquote {
  border: none;
  border-left: 3px solid var(--primary);
  background: transparent;
  padding: 0.1em 0 0.1em 1.2em;
  margin: 1.5em 0;
  color: var(--app-text-secondary);
  font-style: normal;
}

.markdown-body blockquote p:last-child {
  margin-bottom: 0;
}

.markdown-body table {
  border-collapse: separate;
  border-spacing: 0;
  width: 100%;
  margin: 1.5em 0;
  border: 1px solid var(--line);
  border-radius: var(--radius);
  overflow: hidden;
  font-size: 0.94em;
}

.markdown-body th,
.markdown-body td {
  border: none;
  border-bottom: 1px solid var(--line);
  border-right: 1px solid var(--line);
  padding: 0.6em 0.9em;
  text-align: left;
}

.markdown-body th:last-child,
.markdown-body td:last-child {
  border-right: none;
}

.markdown-body tr:last-child td {
  border-bottom: none;
}

.markdown-body tbody tr:hover td {
  background: var(--surface-2);
}

.markdown-body th {
  background-color: var(--app-markdown-code-header-bg);
  font-weight: 600;
  color: var(--app-text-emphasis);
}

.markdown-body img {
  max-width: 100%;
  height: auto;
  display: block;
  margin: 1.5em auto;
  border-radius: var(--radius);
  border: 1px solid var(--line);
}
</style>
