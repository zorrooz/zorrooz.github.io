import { unified } from 'unified'
import remarkParse from 'remark-parse'
import remarkGfm from 'remark-gfm'     // 支持表格、任务列表等
import remarkMath from 'remark-math'
import remarkRehype from 'remark-rehype'
import remarkBreaks from 'remark-breaks' // 新增：处理换行
import rehypeKatex from 'rehype-katex'
import rehypeStringify from 'rehype-stringify'

const processor = unified()
  .use(remarkParse)
  .use(remarkGfm)
  .use(remarkBreaks)        // 新增：启用换行处理
  .use(remarkMath)
  .use(remarkRehype)
  .use(rehypeKatex, {
    throwOnError: false,
    errorColor: '#cc0000'
  })
  .use(rehypeStringify)

export async function renderMarkdown(markdown) {
  const result = await processor.process(markdown)
  return String(result)
}
