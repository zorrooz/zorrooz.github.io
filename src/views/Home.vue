<!-- Home.vue -->
<template>
  <div class="container py-4">
    <div class="row">

      <!-- 中间主内容区域 -->
      <div class="col-12 col-lg-9 order-1 order-lg-2">
        <div class="container">
          <div class="row">
            <div class="col">
              <MD :content="markdownContent" />
              <p>abcd 测试</p>
            </div>
          </div>
        </div>
      </div>

      <!-- 左侧边栏：移动端显示在最下面 -->
      <div class="col-12 col-lg-3 order-2 order-lg-1">
        <div class="d-flex flex-column align-items-center align-items-lg-end w-100 gap-3">
          <ProfileCard class="w-100" />
          <TagCloud class="w-100" />
        </div>
      </div>

    </div>
  </div>
</template>

<script>
import ProfileCard from '@/components/layout/ProfileCard.vue'
import TagCloud from '@/components/layout/TagCloud.vue'
import MD from '@/components/layout/MD.vue'
import markdownContent from '../content/intro.md?raw'

export default {
  name: 'MainView',
  components: {
    ProfileCard,
    TagCloud,
    MD
  },
  data() {
    return {
      markdownContent
    }
  }
}
</script>

<style scoped></style>
