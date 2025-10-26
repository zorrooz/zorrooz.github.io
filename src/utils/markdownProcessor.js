//@/utils/markdownProcessor.js
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
import julia from 'highlight.js/lib/languages/julia'
import dockerfile from 'highlight.js/lib/languages/dockerfile'
import 'highlight.js/styles/github.css'
const languages = { ...common, julia, dockerfile }

const processor = unified()
  .use(remarkParse)
  .use(remarkGfm)
  .use(remarkBreaks)
  .use(remarkMath)
  .use(remarkRehype)
  .use(rehypeHighlight, { languages })
  .use(rehypeKatex, { throwOnError: false, errorColor: '#cc0000' })
  .use(rehypeStringify)
  .use(remarkFrontmatter, ['yaml', 'toml'])
  .use(remarkParseFrontmatter)
export async function renderMarkdown(markdown) {

  const result = await processor.process(markdown)
  return String(result)
}
