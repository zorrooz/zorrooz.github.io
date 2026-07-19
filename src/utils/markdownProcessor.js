import { unified } from 'unified'
import remarkParse from 'remark-parse'
import remarkFrontmatter from 'remark-frontmatter'
import remarkParseFrontmatter from 'remark-parse-frontmatter'
import remarkGfm from 'remark-gfm'
import remarkMath from 'remark-math'
import remarkBreaks from 'remark-breaks'
import remarkRehype from 'remark-rehype'
import rehypeHighlight from 'rehype-highlight'
import rehypeKatex from 'rehype-katex'
import rehypeStringify from 'rehype-stringify'

import { common } from 'lowlight'

const loadHighlightStyle = () => {
  if (typeof window === 'undefined') return

  const isDark = document.documentElement.getAttribute('data-bs-theme') === 'dark'
  const existingStyle = document.getElementById('highlight-style')

  if (existingStyle) {
    existingStyle.remove()
  }

  const link = document.createElement('link')
  link.id = 'highlight-style'
  link.rel = 'stylesheet'
  link.href = isDark
    ? 'https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.11.1/styles/github-dark-dimmed.min.css'
    : 'https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.11.1/styles/github.min.css'

  document.head.appendChild(link)
}

if (typeof window !== 'undefined') {
  loadHighlightStyle()

  const observer = new MutationObserver((mutations) => {
    mutations.forEach((mutation) => {
      if (mutation.attributeName === 'data-bs-theme') {
        loadHighlightStyle()
      }
    })
  })

  observer.observe(document.documentElement, {
    attributes: true,
    attributeFilter: ['data-bs-theme'],
  })
}
const processor = unified()
  .use(remarkParse)
  .use(remarkFrontmatter, ['yaml', 'toml'])
  .use(remarkParseFrontmatter)
  .use(remarkGfm)
  .use(remarkBreaks)
  .use(remarkMath)
  .use(remarkRehype)
  .use(rehypeHighlight, { languages: common })
  .use(rehypeKatex, { throwOnError: false, errorColor: '#cc0000' })
  .use(rehypeStringify)
export async function renderMarkdown(markdown) {
  const result = await processor.process(markdown)
  return String(result)
}
