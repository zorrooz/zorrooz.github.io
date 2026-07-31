// 构建期专用（Node 环境），勿在客户端 import
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

export async function renderMarkdown(markdown: string): Promise<string> {
  const result = await processor.process(markdown)
  return String(result)
}
