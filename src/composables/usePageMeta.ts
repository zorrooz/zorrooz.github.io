/**
 * 统一页面 SEO title：`页面标题 - zorrooz’s blog`。
 * 用法：列表页 `usePageMeta(t('metaHome'))`；文章页 `usePageMeta(computed(() => post.title))`；
 * 不传则只用品牌名。
 */
import { computed, type Ref } from 'vue'
import { useHead } from '@unhead/vue'

import { BLOG_TITLE } from '@/config'

export function usePageMeta(title?: string | Ref<string | undefined>) {
  useHead({
    title: computed(() => {
      const t = typeof title === 'string' ? title : title?.value
      return t ? `${t} - ${BLOG_TITLE}` : BLOG_TITLE
    }),
  })
}
