/**
 * 构建期 markdown → HTML 渲染管线（构建期专用，勿在客户端 import）。
 */
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

/**
 * 链接协议白名单：仅允许 http/https/mailto；相对路径与锚点（无协议）保留。
 * 其余协议（javascript:、data:、vbscript:、tel: 等）移除 href——设计决策。
 */
const SAFE_LINK_PROTOCOL = /^(https?:|mailto:)/i
/**
 * 清理首尾 C0 控制字符与空白：浏览器 URL 解析会剥离它们，须先剥离再检测协议。
 */
// eslint-disable-next-line no-control-regex -- 刻意匹配 C0 控制字符以封堵协议注入绕过
const TRIM_CONTROLS = new RegExp('^[\\u0000-\\u0020]+|[\\u0000-\\u0020]+$', 'g')

function sanitizeLinkHrefs() {
  return (tree: any) => {
    const walk = (node: any) => {
      if (
        node?.type === 'element' &&
        node.tagName === 'a' &&
        typeof node.properties?.href === 'string'
      ) {
        const href = node.properties.href.replace(TRIM_CONTROLS, '')
        const scheme = href.match(/^([a-zA-Z][a-zA-Z0-9+.-]*):/)
        if (scheme && !SAFE_LINK_PROTOCOL.test(href)) delete node.properties.href
      }
      node?.children?.forEach(walk)
    }
    walk(tree)
  }
}

const processor = unified()
  .use(remarkParse)
  .use(remarkFrontmatter, ['yaml', 'toml'])
  .use(remarkParseFrontmatter)
  .use(remarkGfm)
  .use(remarkBreaks)
  .use(remarkMath)
  .use(remarkRehype)
  .use(sanitizeLinkHrefs)
  .use(rehypeHighlight, { languages: common })
  .use(rehypeKatex, { throwOnError: false, errorColor: '#cc0000' })
  .use(rehypeStringify)

export async function renderMarkdown(markdown: string): Promise<string> {
  const result = await processor.process(markdown)
  return String(result)
}
